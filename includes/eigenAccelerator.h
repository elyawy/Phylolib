#ifndef ___EIGEN_ACCELERATOR
#define ___EIGEN_ACCELERATOR

#include "pijAccelerator.h"
#include "replacementModel.h"
#include "fromQtoPt.h"
#include "errorMsg.h"
#include <cmath>
#include <limits>
#include <cassert>

// Accelerates P(t) computation using eigen decomposition of Q.
//
// At construction: decomposes Q into L, eigenvalues λ, R such that
//   P(t) = L · diag(exp(λ*t)) · R
//
// Per distinct t:
//   1. Scale rows of R by exp(λk * t)  (N exp calls)
//   2. P = L · R_scaled               (N×N matrix multiply)
//
// L and R are stored as flat stack-friendly N×N arrays (no heap).
// The full matrix is cached; Pij_t(i,j,t) refills only when t changes.
// N is a compile-time constant, enabling loop unrolling and auto-vectorization.

template<int N>
class eigenAccelerator : public pijAccelerator {
public:
    explicit eigenAccelerator(const replacementModel* pb)
        : _pb(pb->clone()),
          _cachedT(std::numeric_limits<MDOUBLE>::quiet_NaN())
    {
        assert(pb->alphabetSize() == N);
        buildDecomp();
    }

    eigenAccelerator(const eigenAccelerator& other)
        : _pb(other._pb ? other._pb->clone() : nullptr),
          _cachedT(other._cachedT)
    {
        for (int i = 0; i < N; ++i) {
            _eigenVec[i] = other._eigenVec[i];
            for (int j = 0; j < N; ++j) {
                _L[i][j]   = other._L[i][j];
                _R[i][j]   = other._R[i][j];
                _mat[i][j] = other._mat[i][j];
            }
        }
    }

    const MDOUBLE Qij(const int i, const int j) const { return _pb->Qij(i, j); }

    inline const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const {
        if (t != _cachedT) refill(t);
        return _mat[i][j];
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

    void refill(const MDOUBLE t) const
    {
        // Step 1: R_scaled[k][j] = R[k][j] * exp(λk * t)
        double R_scaled[N][N];
        for (int k = 0; k < N; ++k) {
            const double ek = std::exp(_eigenVec[k] * t);
            for (int j = 0; j < N; ++j)
                R_scaled[k][j] = _R[k][j] * ek;
        }

        // Step 2: P = L * R_scaled
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                double sum = 0.0;
                for (int k = 0; k < N; ++k)
                    sum += _L[i][k] * R_scaled[k][j];
                _mat[i][j] = sum;
            }

        _cachedT = t;
    }

    replacementModel* _pb;

    // Eigen decomposition — flat 2D arrays, stack-friendly
    double _eigenVec[N];     // eigenvalues λ[k]
    double _L[N][N];         // left eigenvectors  L[i][k]
    double _R[N][N];         // right eigenvectors R[k][j]

    // Cached P(t) matrix and the t it was computed for
    mutable double  _mat[N][N];
    mutable MDOUBLE _cachedT;
};

#endif