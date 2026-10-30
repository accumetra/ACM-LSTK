/*=========================================================================

  Sub-voxel edge refinement of an LSTK segmentation (--RefineEdge).

  LSTK's level set is used only as a guide: where the lesion is, and where it is cut from attached structures.

  1. Intensities. The lung level is the median of the CT in a shell 2-4 mm outside LSTK's surface. The lesion level is
     the median of the core (>= 1 mm, or 80% of the deepest point, inside) of the half-max region grown from the seed.
     The threshold T is halfway between them: the half maximum of the edge profile.
  2. FWHM surface. On the level set's grid, the CT (lightly smoothed, cubic B-spline interpolated) is thresholded at T,
     keeping the part connected to the seed. If that part runs out of the working crop (Margin mm around LSTK's
     region), something large and dense is attached (vessel, wall): growth is then limited to LSTK's region plus Band
     mm, and not within Free mm of the dense structure beyond it, so LSTK's cut is kept there. Otherwise the lesion is
     isolated and the whole part is used. The marching-cubes surface of the thresholded CT is the sub-voxel FWHM
     surface.
  3. Blur and curvature. Blurring moves the half-max point of a curved edge inward by 1/2 tr(K P S P), where K is the
     curvature of the CT's level surface through the point (Gaussian derivatives), P = I - n n^T and S the blur
     covariance. S is the scanner PSF if given, otherwise it is fitted from the edge slopes: at the half-max point the
     slope is (lesion - lung) / (sqrt(2 pi) sigma_n), with sigma_n^2 = n^T S n. Each surface point is treated as a
     sphere of the same mean curvature, and the exact half-max radius of a blurred sphere is inverted to give its
     outward shift. The corrected volume is the FWHM volume plus the sum of shift x area.

  Tested on simulated ellipsoids (ACM-Simulator, random orientations and sub-voxel positions): within about +1% of
  the true volume from 8 mm diameters up, against -10% to +900% for the level set alone.

=========================================================================*/
#ifndef SubvoxelEdgeRefinement_h
#define SubvoxelEdgeRefinement_h

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"

namespace lstk
{

class SubvoxelEdgeRefinement
{
public:
  using CTImageType = itk::Image< short, 3 >;
  using RealImageType = itk::Image< float, 3 >;
  using MaskImageType = itk::Image< unsigned char, 3 >;
  using PointType = RealImageType::PointType;
  using IndexType = RealImageType::IndexType;
  using RegionType = RealImageType::RegionType;

  // Parameters (mm)
  double Band = 0.6;        // growth allowed beyond LSTK's surface next to an attachment
  double Smooth = 0.3;      // Gaussian applied to the CT before interpolation (counted in the blur)
  double Margin = 8.0;      // working crop around LSTK's region; a half-max region reaching its edge is attached
  double Free = 1.0;        // no growth, and no curvature shift, this close to an attached structure
  double CurvatureSigma = 0.4;
  double LevelSetIso = -0.5;  // LSTK's surface in its level set (its region is level set > LevelSetIso)

  // Results
  bool Valid = false;
  bool Attached = false;
  double Lesion = 0, Lung = 0, Threshold = 0;
  double FWHMVolume = 0, CorrectedVolume = 0;
  double SigmaUsed[3] = { 0, 0, 0 };          // blur used for the correction (x, y, z), including Smooth
  double SigmaEstimated[3] = { 0, 0, 0 };     // blur fitted from the edges (x, y, z), including Smooth
  RealImageType::Pointer FWHMLevelSet;        // same grid as the input level set; surface at -0.5
  RealImageType::Pointer CorrectedLevelSet;

  // psfSigma: scanner PSF sigma (x, y, z) in mm, or nullptr to estimate it from the edges
  void Run(const CTImageType * ct, const RealImageType * levelSet, const PointType & seed, const double * psfSigma)
  {
    m_LevelSet = levelSet;
    const double h = levelSet->GetSpacing()[0];
    m_H = h;

    // ---- working crop on the level set grid around LSTK's region and the seed
    RegionType full = levelSet->GetLargestPossibleRegion();
    IndexType lo, hi;
    for (int d = 0; d < 3; d++) { lo[d] = std::numeric_limits< long >::max(); hi[d] = std::numeric_limits< long >::min(); }
    bool any = false;
    for (itk::ImageRegionConstIteratorWithIndex< RealImageType > it(levelSet, full); !it.IsAtEnd(); ++it)
      if (it.Get() > LevelSetIso)
      {
        any = true;
        for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], it.GetIndex()[d]); hi[d] = std::max(hi[d], it.GetIndex()[d]); }
      }
    IndexType seedIndex;
    levelSet->TransformPhysicalPointToIndex(seed, seedIndex);
    if (!any) { lo = seedIndex; hi = seedIndex; }
    const long pad = static_cast< long >(std::ceil(m_Margin() / h));
    IndexType cropStart; RealImageType::SizeType cropSize;
    for (int d = 0; d < 3; d++)
    {
      const long a = std::max< long >(std::min(lo[d], seedIndex[d]) - pad, full.GetIndex()[d]);
      const long b = std::min< long >(std::max(hi[d], seedIndex[d]) + pad, full.GetIndex()[d] + full.GetSize()[d] - 1);
      cropStart[d] = a; cropSize[d] = b - a + 1;
    }
    m_Crop = RegionType(cropStart, cropSize);

    // ---- CT: extract around the crop, smooth, cubic B-spline onto the crop grid
    PointType p0, p1;
    levelSet->TransformIndexToPhysicalPoint(cropStart, p0);
    IndexType cropEnd; for (int d = 0; d < 3; d++) cropEnd[d] = cropStart[d] + cropSize[d] - 1;
    levelSet->TransformIndexToPhysicalPoint(cropEnd, p1);
    CTImageType::IndexType c0, c1;
    ct->TransformPhysicalPointToIndex(p0, c0); ct->TransformPhysicalPointToIndex(p1, c1);
    CTImageType::RegionType ctRegion;
    for (int d = 0; d < 3; d++)
    {
      const long m = 4 + static_cast< long >(std::ceil(4 * Smooth / ct->GetSpacing()[d]));
      const long a = std::min(c0[d], c1[d]) - m, b = std::max(c0[d], c1[d]) + m;
      ctRegion.SetIndex(d, a); ctRegion.SetSize(d, b - a + 1);
    }
    ctRegion.Crop(ct->GetLargestPossibleRegion());
    using ROIType = itk::RegionOfInterestImageFilter< CTImageType, CTImageType >;
    auto roi = ROIType::New(); roi->SetInput(ct); roi->SetRegionOfInterest(ctRegion);
    using SmoothType = itk::DiscreteGaussianImageFilter< CTImageType, RealImageType >;
    auto smooth = SmoothType::New(); smooth->SetInput(roi->GetOutput());
    smooth->SetVariance(Smooth * Smooth); smooth->SetUseImageSpacing(true); smooth->SetMaximumKernelWidth(64);
    smooth->Update();
    using SplineType = itk::BSplineInterpolateImageFunction< RealImageType, double, double >;
    auto spline = SplineType::New(); spline->SetSplineOrder(3); spline->SetInputImage(smooth->GetOutput());

    m_I = NewReal();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(m_I, m_Crop); !it.IsAtEnd(); ++it)
    {
      PointType p; m_I->TransformIndexToPhysicalPoint(it.GetIndex(), p);
      itk::ContinuousIndex< double, 3 > ci;
      smooth->GetOutput()->TransformPhysicalPointToContinuousIndex(p, ci);
      const auto & r = smooth->GetOutput()->GetLargestPossibleRegion();
      for (int d = 0; d < 3; d++)
        ci[d] = std::min< double >(std::max< double >(ci[d], r.GetIndex()[d]), r.GetIndex()[d] + r.GetSize()[d] - 1);
      it.Set(static_cast< float >(spline->EvaluateAtContinuousIndex(ci)));
    }

    // ---- 1. intensities
    auto M = NewMask();
    for (itk::ImageRegionIteratorWithIndex< MaskImageType > it(M, m_Crop); !it.IsAtEnd(); ++it)
      it.Set(levelSet->GetPixel(it.GetIndex()) > LevelSetIso ? 1 : 0);
    auto dout = OutsideDistance(M);                    // distance to LSTK's region (0 inside)
    std::vector< float > shell;
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(dout, m_Crop); !it.IsAtEnd(); ++it)
      if (it.Get() >= 2 && it.Get() <= 4) shell.push_back(m_I->GetPixel(it.GetIndex()));
    if (shell.empty()) return;
    Lung = Median(shell);

    auto D = NewMask();
    for (itk::ImageRegionIteratorWithIndex< MaskImageType > it(D, m_Crop); !it.IsAtEnd(); ++it)
      it.Set(dout->GetPixel(it.GetIndex()) <= Band ? 1 : 0);
    m_Seed = seedIndex;
    if (!m_Crop.IsInside(m_Seed)) return;

    std::vector< float > ball;
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(m_I, m_Crop); !it.IsAtEnd(); ++it)
    {
      PointType p; m_I->TransformIndexToPhysicalPoint(it.GetIndex(), p);
      if (p.SquaredEuclideanDistanceTo(seed) <= 1.0) ball.push_back(it.Get());
    }
    if (ball.empty()) return;
    Lesion = Median(ball);
    for (int iter = 0; iter < 2; iter++)
    {
      auto R = Region(0.5 * (Lesion + Lung), D);
      auto din = InsideDistance(R);
      float dmax = 0;
      for (itk::ImageRegionIterator< RealImageType > it(din, m_Crop); !it.IsAtEnd(); ++it) dmax = std::max(dmax, it.Get());
      const float depth = std::min(1.0f, 0.8f * dmax);
      std::vector< float > core;
      for (itk::ImageRegionIteratorWithIndex< RealImageType > it(din, m_Crop); !it.IsAtEnd(); ++it)
        if (it.Get() >= depth && it.Get() > 0) core.push_back(m_I->GetPixel(it.GetIndex()));
      if (core.size() >= 10) Lesion = Median(core);
    }
    Threshold = 0.5 * (Lesion + Lung);

    auto all = NewMask(); all->FillBuffer(1);
    auto E = Region(Threshold, all);
    Attached = TouchesCropEdge(E);
    auto allowed = NewMask(); allowed->FillBuffer(1);
    RealImageType::Pointer dX;
    MaskImageType::Pointer R;
    if (Attached)
    {
      auto X = NewMask();     // dense structure beyond the band
      for (itk::ImageRegionIteratorWithIndex< MaskImageType > it(X, m_Crop); !it.IsAtEnd(); ++it)
        it.Set(E->GetPixel(it.GetIndex()) && !D->GetPixel(it.GetIndex()) ? 1 : 0);
      dX = OutsideDistance(X);
      for (itk::ImageRegionIteratorWithIndex< MaskImageType > it(allowed, m_Crop); !it.IsAtEnd(); ++it)
      {
        const auto & i = it.GetIndex();
        it.Set(M->GetPixel(i) || (D->GetPixel(i) && dX->GetPixel(i) > Free) ? 1 : 0);
      }
      R = Region(Threshold, allowed);
    }
    else
      R = E;

    // ---- 2. FWHM surface: threshold function, forced below T away from the seed's component
    auto F = NewReal();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(F, m_Crop); !it.IsAtEnd(); ++it)
    {
      const auto & i = it.GetIndex();
      const float f = m_I->GetPixel(i) - static_cast< float >(Threshold);
      it.Set(DilatedAt(R, i) && allowed->GetPixel(i) ? f : std::min(f, -1.0f));
    }
    std::vector< std::array< double, 3 > > V;
    std::vector< std::array< vtkIdType, 3 > > T;
    MarchingCubes(F, V, T);
    if (T.empty()) return;
    std::vector< double > area(V.size(), 0.0);
    FWHMVolume = VolumeAndVertexArea(V, T, area);

    // ---- 3. blur estimate from edge slopes, and curvature shift
    const size_t nv = V.size();
    std::vector< std::array< double, 3 > > n(nv);
    std::vector< double > gm(nv);
    for (size_t k = 0; k < nv; k++)
    {
      std::array< double, 3 > g = Trilinear3(V[k], [this](const IndexType & i) { return CentralGradient(i); });
      gm[k] = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
      for (int d = 0; d < 3; d++) n[k][d] = -g[d] / std::max(gm[k], 1e-6);
    }
    std::vector< double > sorted(gm); std::nth_element(sorted.begin(), sorted.begin() + nv / 2, sorted.end());
    const double gmed = sorted[nv / 2];
    std::vector< char > ok(nv);
    for (size_t k = 0; k < nv; k++)
    {
      ok[k] = gm[k] > 0.5 * gmed;
      if (Attached && ok[k])
        ok[k] = Trilinear3(V[k], [&dX](const IndexType & i) {
                  const double v = std::min(99.0, static_cast< double >(dX->GetPixel(i))); return std::array< double, 3 >{ v, 0, 0 }; })[0]
                > Free + 0.5;
    }
    // weighted least squares: sigma_n^2 = sum_d n_d^2 S_d
    double A[3][3] = { { 0 } }, b[3] = { 0 };
    for (size_t k = 0; k < nv; k++)
    {
      if (!ok[k]) continue;
      const double s2 = std::pow((Lesion - Lung) / (std::sqrt(2 * M_PI) * std::max(gm[k], 1e-6)), 2);
      const double a[3] = { n[k][0] * n[k][0], n[k][1] * n[k][1], n[k][2] * n[k][2] };
      for (int i = 0; i < 3; i++) { b[i] += area[k] * a[i] * s2; for (int j = 0; j < 3; j++) A[i][j] += area[k] * a[i] * a[j]; }
    }
    double Sest[3];
    Solve3(A, b, Sest);
    for (int d = 0; d < 3; d++) { Sest[d] = std::max(0.0, Sest[d]); SigmaEstimated[d] = std::sqrt(Sest[d]); }
    double S[3];
    for (int d = 0; d < 3; d++) S[d] = psfSigma ? psfSigma[d] * psfSigma[d] + Smooth * Smooth : Sest[d];
    for (int d = 0; d < 3; d++) SigmaUsed[d] = std::sqrt(S[d]);

    // level-surface curvature from Gaussian derivatives of the CT, at the vertices and in a band of voxels
    std::vector< IndexType > bandVoxels;
    std::vector< double > bandDist;
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(F, m_Crop); !it.IsAtEnd(); ++it)
    {
      const auto g = CentralGradient(it.GetIndex());
      const double gn = std::max(std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]), 0.25 * gmed);
      const double dist = it.Get() / gn;    // approximate signed distance to the FWHM surface, positive inside
      if (std::fabs(dist) < 1.5) { bandVoxels.push_back(it.GetIndex()); bandDist.push_back(dist); }
    }
    std::vector< std::array< double, 9 > > dv(nv), db(bandVoxels.size());   // gx gy gz hxx hyy hzz hxy hxz hyz
    const int orders[9][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 },
                               { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 } };
    for (int q = 0; q < 9; q++)
    {
      auto Dq = GaussianDerivative(m_I, orders[q], CurvatureSigma);
      for (size_t k = 0; k < nv; k++)
        dv[k][q] = Trilinear3(V[k], [&Dq](const IndexType & i) { return std::array< double, 3 >{ Dq->GetPixel(i), 0, 0 }; })[0];
      for (size_t k = 0; k < bandVoxels.size(); k++) db[k][q] = Dq->GetPixel(bandVoxels[k]);
    }

    double shiftSum = 0;
    for (size_t k = 0; k < nv; k++)
      if (ok[k]) shiftSum += Shift(dv[k], n[k], S) * area[k];
    CorrectedVolume = FWHMVolume + shiftSum;

    // ---- output level sets on the input grid (surface at -0.5, like LSTK's)
    // Dense voxels of other structures (above T, outside the seed's region or the allowed growth) are set well outside;
    // voxels below T keep their distance, so the shift can move the surface into them
    const auto kept = [&](const IndexType & i) {
      return allowed->GetPixel(i) && (DilatedAt(R, i) || m_I->GetPixel(i) <= Threshold); };
    FWHMLevelSet = NewFullReal(); CorrectedLevelSet = NewFullReal();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(F, m_Crop); !it.IsAtEnd(); ++it)
    {
      const auto g = CentralGradient(it.GetIndex());
      const double gn = std::max(std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]), 0.25 * gmed);
      const double dist = kept(it.GetIndex()) ? std::max(-4.0, std::min(4.0, it.Get() / gn)) : -4.0;
      FWHMLevelSet->SetPixel(it.GetIndex(), static_cast< float >(dist - 0.5));
      CorrectedLevelSet->SetPixel(it.GetIndex(), static_cast< float >(dist - 0.5));
    }
    for (size_t k = 0; k < bandVoxels.size(); k++)
    {
      const auto & i = bandVoxels[k];
      if (!kept(i)) continue;
      const auto g = CentralGradient(i);
      const double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
      if (gn <= 0.5 * gmed) continue;
      if (Attached && dX->GetPixel(i) <= Free + 0.5) continue;
      std::array< double, 3 > nn = { -g[0] / gn, -g[1] / gn, -g[2] / gn };
      CorrectedLevelSet->SetPixel(i, static_cast< float >(bandDist[k] + Shift(db[k], nn, S) - 0.5));
    }
    KeepSeedComponent(FWHMLevelSet); KeepSeedComponent(CorrectedLevelSet);
    Valid = true;
  }

  // Volume inside the surface at -0.5 of a level set, from marching cubes (for checking the output images)
  static double LevelSetVolume(const RealImageType * ls)
  {
    SubvoxelEdgeRefinement r; r.m_LevelSet = ls; r.m_H = ls->GetSpacing()[0]; r.m_Crop = ls->GetLargestPossibleRegion();
    auto F = r.NewReal();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(F, r.m_Crop); !it.IsAtEnd(); ++it) it.Set(ls->GetPixel(it.GetIndex()) + 0.5f);
    std::vector< std::array< double, 3 > > V; std::vector< std::array< vtkIdType, 3 > > T;
    r.MarchingCubes(F, V, T);
    std::vector< double > area(V.size());
    return T.empty() ? 0.0 : VolumeAndVertexArea(V, T, area);
  }

private:
  const RealImageType * m_LevelSet = nullptr;
  RegionType m_Crop;
  IndexType m_Seed;
  double m_H = 0.2;
  RealImageType::Pointer m_I;

  double m_Margin() const { return Margin; }

  RealImageType::Pointer NewReal() const
  {
    auto im = RealImageType::New(); im->CopyInformation(m_LevelSet); im->SetRegions(m_Crop); im->Allocate(); im->FillBuffer(0); return im;
  }
  RealImageType::Pointer NewFullReal() const
  {
    auto im = RealImageType::New(); im->CopyInformation(m_LevelSet); im->SetRegions(m_LevelSet->GetLargestPossibleRegion());
    im->Allocate(); im->FillBuffer(-4.0f); return im;
  }
  MaskImageType::Pointer NewMask() const
  {
    auto im = MaskImageType::New(); im->CopyInformation(m_LevelSet); im->SetRegions(m_Crop); im->Allocate(); im->FillBuffer(0); return im;
  }

  static float Median(std::vector< float > & v)
  {
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
  }

  // Distance (mm) from each voxel to the nearest voxel of the mask; 0 inside it. Very large if the mask is empty.
  RealImageType::Pointer OutsideDistance(MaskImageType * m) const
  {
    auto out = NewReal();
    bool any = false;
    for (itk::ImageRegionIterator< MaskImageType > it(m, m_Crop); !it.IsAtEnd() && !any; ++it) any = it.Get();
    if (!any) { out->FillBuffer(1e9f); return out; }
    using DistType = itk::SignedMaurerDistanceMapImageFilter< MaskImageType, RealImageType >;
    auto dist = DistType::New(); dist->SetInput(m); dist->SetUseImageSpacing(true); dist->SetSquaredDistance(false);
    dist->SetInsideIsPositive(false); dist->SetBackgroundValue(0); dist->Update();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(out, m_Crop); !it.IsAtEnd(); ++it)
      it.Set(m->GetPixel(it.GetIndex()) ? 0.0f : dist->GetOutput()->GetPixel(it.GetIndex()));
    return out;
  }

  // Distance (mm) from each voxel of the mask to the nearest voxel outside it; 0 outside.
  RealImageType::Pointer InsideDistance(MaskImageType * m) const
  {
    auto inv = NewMask();
    for (itk::ImageRegionIteratorWithIndex< MaskImageType > it(inv, m_Crop); !it.IsAtEnd(); ++it) it.Set(m->GetPixel(it.GetIndex()) ? 0 : 1);
    return OutsideDistance(inv);
  }

  // Only the inside (> -0.5) part face-connected to the seed stays inside; the rest of the crop is set outside
  void KeepSeedComponent(RealImageType * ls) const
  {
    using CTType = itk::ConnectedThresholdImageFilter< RealImageType, MaskImageType >;
    auto ctf = CTType::New(); ctf->SetInput(ls);
    ctf->SetLower(-0.5f + 1e-6f); ctf->SetUpper(std::numeric_limits< float >::max());
    ctf->SetReplaceValue(1); ctf->AddSeed(m_Seed); ctf->Update();
    MaskImageType::Pointer keep = ctf->GetOutput();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(ls, m_Crop); !it.IsAtEnd(); ++it)
      if (it.Get() > -0.5f && !keep->GetPixel(it.GetIndex())) it.Set(-4.5f);
  }

  // Part of (I > T) inside 'allowed' that is face-connected to the seed
  MaskImageType::Pointer Region(double T, MaskImageType * allowed) const
  {
    auto J = NewReal();
    for (itk::ImageRegionIteratorWithIndex< RealImageType > it(J, m_Crop); !it.IsAtEnd(); ++it)
      it.Set(allowed->GetPixel(it.GetIndex()) ? m_I->GetPixel(it.GetIndex()) : -1e30f);
    using CTType = itk::ConnectedThresholdImageFilter< RealImageType, MaskImageType >;
    auto ctf = CTType::New(); ctf->SetInput(J);
    ctf->SetLower(static_cast< float >(T) + 1e-3f); ctf->SetUpper(std::numeric_limits< float >::max());
    ctf->SetReplaceValue(1); ctf->AddSeed(m_Seed); ctf->Update();
    MaskImageType::Pointer r = ctf->GetOutput(); r->DisconnectPipeline();
    return r;
  }

  bool TouchesCropEdge(MaskImageType * m) const
  {
    for (itk::ImageRegionIteratorWithIndex< MaskImageType > it(m, m_Crop); !it.IsAtEnd(); ++it)
    {
      if (!it.Get()) continue;
      for (int d = 0; d < 3; d++)
        if (it.GetIndex()[d] == m_Crop.GetIndex()[d] || it.GetIndex()[d] == m_Crop.GetIndex()[d] + long(m_Crop.GetSize()[d]) - 1) return true;
    }
    return false;
  }

  bool DilatedAt(MaskImageType * m, const IndexType & i) const    // 6-connected dilation by one voxel
  {
    if (m->GetPixel(i)) return true;
    for (int d = 0; d < 3; d++)
      for (int s = -1; s <= 1; s += 2)
      {
        IndexType j = i; j[d] += s;
        if (m_Crop.IsInside(j) && m->GetPixel(j)) return true;
      }
    return false;
  }

  std::array< double, 3 > CentralGradient(const IndexType & i) const   // like numpy.gradient (one-sided at the edge)
  {
    std::array< double, 3 > g;
    for (int d = 0; d < 3; d++)
    {
      IndexType a = i, b = i; a[d] -= 1; b[d] += 1;
      const bool ia = m_Crop.IsInside(a), ib = m_Crop.IsInside(b);
      const double fa = m_I->GetPixel(ia ? a : i), fb = m_I->GetPixel(ib ? b : i);
      g[d] = (fb - fa) / (m_H * ((ia ? 1 : 0) + (ib ? 1 : 0)));
    }
    return g;
  }

  // Trilinear interpolation at a physical point of a per-voxel function returning up to 3 values
  template < class Fn > std::array< double, 3 > Trilinear3(const std::array< double, 3 > & p, Fn f) const
  {
    PointType pt; pt[0] = p[0]; pt[1] = p[1]; pt[2] = p[2];
    itk::ContinuousIndex< double, 3 > ci;
    m_I->TransformPhysicalPointToContinuousIndex(pt, ci);
    IndexType base; double t[3];
    for (int d = 0; d < 3; d++)
    {
      const long lo = m_Crop.GetIndex()[d], hi = lo + long(m_Crop.GetSize()[d]) - 2;
      base[d] = std::min(std::max(static_cast< long >(std::floor(ci[d])), lo), hi);
      t[d] = std::min(1.0, std::max(0.0, ci[d] - base[d]));
    }
    std::array< double, 3 > out = { 0, 0, 0 };
    for (int c = 0; c < 8; c++)
    {
      IndexType j = base; double w = 1;
      for (int d = 0; d < 3; d++) { const int bit = (c >> d) & 1; j[d] += bit; w *= bit ? t[d] : 1 - t[d]; }
      if (w == 0) continue;
      const auto v = f(j);
      for (int q = 0; q < 3; q++) out[q] += w * v[q];
    }
    return out;
  }

  static RealImageType::Pointer GaussianDerivative(RealImageType * in, const int order[3], double sigma)
  {
    using RGType = itk::RecursiveGaussianImageFilter< RealImageType, RealImageType >;
    RealImageType::Pointer cur = in;
    for (unsigned int d = 0; d < 3; d++)
    {
      auto f = RGType::New(); f->SetInput(cur); f->SetDirection(d); f->SetSigma(sigma);
      f->SetOrder(order[d] == 0 ? itk::GaussianOrderEnum::ZeroOrder
                  : order[d] == 1 ? itk::GaussianOrderEnum::FirstOrder : itk::GaussianOrderEnum::SecondOrder);
      f->Update(); cur = f->GetOutput(); cur->DisconnectPipeline();
    }
    // the recursive Gaussian's derivatives are per pixel: convert to per mm
    double scale = 1;
    for (unsigned int d = 0; d < 3; d++) scale *= std::pow(in->GetSpacing()[d], order[d]);
    for (itk::ImageRegionIterator< RealImageType > it(cur, cur->GetBufferedRegion()); !it.IsAtEnd(); ++it)
      it.Set(static_cast< float >(it.Get() / scale));
    return cur;
  }

  // Outward shift from the half-max surface to the true edge at a point with Gaussian derivatives d (gx gy gz hxx hyy
  // hzz hxy hxz hyz) of the CT, outward normal n and blur covariance S (diagonal).
  static double Shift(const std::array< double, 9 > & d, const std::array< double, 3 > & n, const double S[3])
  {
    const double gn = std::max(std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), 1e-6);
    const double H[3][3] = { { d[3], d[6], d[7] }, { d[6], d[4], d[8] }, { d[7], d[8], d[5] } };
    double P[3][3], PSP[3][3];
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) P[a][b] = (a == b) - n[a] * n[b];
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++)
    {
      PSP[a][b] = 0;
      for (int c = 0; c < 3; c++) PSP[a][b] += P[a][c] * S[c] * P[c][b];
    }
    double first = 0, trace = 0;
    for (int a = 0; a < 3; a++) { trace += PSP[a][a]; for (int b = 0; b < 3; b++) first += 0.5 * (-H[a][b] / gn) * PSP[b][a]; }
    return SphereShift(first, trace / 2);
  }

  // Exact half-max radius of a ball of radius R blurred by an isotropic Gaussian of sigma 1
  static double SphereHalfRadius(double R)
  {
    const double a = std::sqrt(2.0);
    auto f = [&](double r) {
      return 0.5 * (std::erf((R - r) / a) + std::erf((R + r) / a)) -
             1.0 / (r * std::sqrt(2 * M_PI)) * (std::exp(-(r - R) * (r - R) / 2) - std::exp(-(r + R) * (r + R) / 2)) - 0.5;
    };
    double lo = 1e-6, hi = R + 6;
    if (f(lo) <= 0) return 0.0;
    for (int i = 0; i < 100; i++) { const double m = 0.5 * (lo + hi); (f(m) > 0 ? lo : hi) = m; }
    return 0.5 * (lo + hi);
  }

  // first = 1/2 tr(K P S P) (first-order shift), sigma2 = tangential variance; the equivalent sphere has half-max
  // radius sigma2 / first. Flat or concave points keep the first-order shift.
  static double SphereShift(double first, double sigma2)
  {
    static std::vector< double > X, Y;
    if (X.empty())
    {
      for (int i = 0; i < 400; i++) Y.push_back(0.8 + i * (10.0 - 0.8) / 399);
      for (int i = 0; i < 400; i++) Y.push_back(10.05 + i * (400.0 - 10.05) / 399);
      for (double y : Y) X.push_back(SphereHalfRadius(y));
    }
    if (first <= 1e-9 || sigma2 <= 0) return first;
    const double s = std::sqrt(sigma2), x = sigma2 / first / s;
    double y;
    if (x > X.back()) y = x + 1 / x;
    else
    {
      const size_t i = std::lower_bound(X.begin(), X.end(), x) - X.begin();
      if (i == 0) y = Y[0];
      else if (X[i] == X[i - 1]) y = Y[i];
      else y = Y[i - 1] + (Y[i] - Y[i - 1]) * (x - X[i - 1]) / (X[i] - X[i - 1]);
    }
    return s * (y - x);
  }

  static void Solve3(const double A[3][3], const double b[3], double x[3])
  {
    auto det = [](const double M[3][3]) {
      return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
             M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    };
    const double d = det(A);
    for (int c = 0; c < 3; c++)
    {
      double Mc[3][3];
      for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Mc[i][j] = (j == c) ? b[i] : A[i][j];
      x[c] = std::fabs(d) > 1e-30 ? det(Mc) / d : 0.0;
    }
  }

  // Marching cubes of F at 0 (physical coordinates; the level set grid has an identity direction)
  void MarchingCubes(RealImageType * F, std::vector< std::array< double, 3 > > & V, std::vector< std::array< vtkIdType, 3 > > & T) const
  {
    auto img = vtkSmartPointer< vtkImageData >::New();
    const auto sz = m_Crop.GetSize();
    img->SetDimensions(sz[0], sz[1], sz[2]);
    PointType o; F->TransformIndexToPhysicalPoint(m_Crop.GetIndex(), o);
    img->SetOrigin(o[0], o[1], o[2]);
    img->SetSpacing(F->GetSpacing()[0], F->GetSpacing()[1], F->GetSpacing()[2]);
    auto arr = vtkSmartPointer< vtkFloatArray >::New();
    arr->SetNumberOfTuples(m_Crop.GetNumberOfPixels());
    vtkIdType k = 0;
    for (itk::ImageRegionIterator< RealImageType > it(F, m_Crop); !it.IsAtEnd(); ++it) arr->SetValue(k++, it.Get());   // x fastest, as VTK
    img->GetPointData()->SetScalars(arr);
    auto mc = vtkSmartPointer< vtkMarchingCubes >::New();
    mc->SetInputData(img); mc->SetValue(0, 0.0); mc->ComputeNormalsOff(); mc->ComputeGradientsOff(); mc->ComputeScalarsOff();
    mc->Update();
    vtkPolyData * pd = mc->GetOutput();
    V.resize(pd->GetNumberOfPoints());
    for (vtkIdType i = 0; i < pd->GetNumberOfPoints(); i++) { double p[3]; pd->GetPoint(i, p); V[i] = { p[0], p[1], p[2] }; }
    auto ids = vtkSmartPointer< vtkIdList >::New();
    vtkCellArray * polys = pd->GetPolys();
    polys->InitTraversal();
    while (polys->GetNextCell(ids))
      if (ids->GetNumberOfIds() == 3) T.push_back({ ids->GetId(0), ids->GetId(1), ids->GetId(2) });
  }

  static double VolumeAndVertexArea(const std::vector< std::array< double, 3 > > & V,
                                    const std::vector< std::array< vtkIdType, 3 > > & T, std::vector< double > & area)
  {
    double vol = 0;
    for (const auto & t : T)
    {
      const auto & a = V[t[0]]; const auto & b = V[t[1]]; const auto & c = V[t[2]];
      vol += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
      const double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, w[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
      const double cx = u[1] * w[2] - u[2] * w[1], cy = u[2] * w[0] - u[0] * w[2], cz = u[0] * w[1] - u[1] * w[0];
      const double A = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
      for (int q = 0; q < 3; q++) area[t[q]] += A / 3;
    }
    return std::fabs(vol) / 6;
  }
};

} // namespace lstk

#endif
