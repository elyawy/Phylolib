#ifndef ___EIGEN_ACCELERATOR
#define ___EIGEN_ACCELERATOR

#include "pijAccelerator.h"
#include "replacementModel.h"
#include "fromQtoPt.h"
#include "errorMsg.h"
#include <cmath>
#include <limits>
#include <cassert>

template<int N>
class eigenAccelerator : public pijAccelerator {
public:
    explicit eigenAccelerator(const replacementModel* pb)
        : _pb(pb->clone()),
          _cachedT(std::numeric_limits<MDOUBLE>::quiet_NaN()),
          _cachedRowIdx(-1)
    {
        assert(pb->alphabetSize() == N);
        buildDecomp();
    }

    eigenAccelerator(const eigenAccelerator& other)
        : _pb(other._pb ? other._pb->clone() : nullptr),
          _cachedT(other._cachedT),
          _cachedRowIdx(other._cachedRowIdx)
    {
        for (int i = 0; i < N; ++i) {
            _eigenVec[i] = other._eigenVec[i];
            _cachedExp[i] = other._cachedExp[i];
            _cachedRow[i] = other._cachedRow[i];
            for (int j = 0; j < N; ++j) {
                _L[i][j]   = other._L[i][j];
                _R[i][j]   = other._R[i][j];
            }
        }
    }

    const MDOUBLE Qij(const int i, const int j) const { return _pb->Qij(i, j); }

    void fillPijRow(const int i, const double t, double* outRow) const override {
        // 1. Refresh exponent cache if the branch time/rate state changes
        if (t != _cachedT) {
            for (int k = 0; k < N; ++k) {
                _cachedExp[k] = std::exp(_eigenVec[k] * t);
            }
            _cachedT = t;
        }

        // 2. Clear outRow buffer
        for (int col = 0; col < N; ++col) {
            outRow[col] = 0.0;
        }

        // 3. Vectorized row math: O(N^2) instead of O(N^3)
        const double* const row_L = _L[i];
        for (int k = 0; k < N; ++k) {
            const double L_scaled = row_L[k] * _cachedExp[k];
            const double* const row_R = _R[k];
            
            // Loop unrolls and vectorizes flawlessly with SIMD alignment
            for (int col = 0; col < N; ++col) {
                outRow[col] += L_scaled * row_R[col];
            }
        }
    }

    // Preserve scalar compliance for backward compatibility 
    inline const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const override {
        double tmpRow[N];
        fillPijRow(i, t, tmpRow);
        return tmpRow[j];
    }
    const MDOUBLE dPij_dt  (const int i, const int j, const MDOUBLE t) const { return _pb->dPij_dt(i, j, t); }
    const MDOUBLE d2Pij_dt2(const int i, const int j, const MDOUBLE t) const { return _pb->d2Pij_dt2(i, j, t); }

    const MDOUBLE freq(const size_t i) const { return _pb->freq(i); }
    virtual pijAccelerator* clone() const { return new eigenAccelerator(*this); }
    virtual ~eigenAccelerator() { delete _pb; }
    virtual replacementModel* getReplacementModel() const { return _pb; }
    virtual const size_t alphabetSize() const { return N; }

private:
    void buildDecomp()
    {
        Vdouble freq(N);
        for (int i = 0; i < N; ++i) freq[i] = _pb->freq(i);

        VVdouble Q(N, Vdouble(N));
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                Q[i][j] = _pb->Qij(i, j);

        q2pt q;
        q.fillFromRateMatrix(freq, Q);

        const VVdouble& L  = q.getLeftEigen();
        const VVdouble& R  = q.getRightEigen();
        const Vdouble&  ev = q.getEigenVec();

        for (int k = 0; k < N; ++k) {
            _eigenVec[k] = ev[k];
            for (int i = 0; i < N; ++i) _L[i][k] = L[i][k];
            for (int j = 0; j < N; ++j) _R[k][j] = R[k][j];
        }
    }

    replacementModel* _pb;
    double _eigenVec[N];
    double _L[N][N];
    double _R[N][N];

    mutable double  _cachedExp[N];
    mutable MDOUBLE _cachedT;

    // Row buffer states
    mutable double  _cachedRow[N];
    mutable int     _cachedRowIdx;
};

#endif