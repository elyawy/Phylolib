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
        // 1D local array guarantees no pointer aliasing
        double tmp[N * N] = {0.0};

        // Outer product loop order: k, i, j
        for (int k = 0; k < N; ++k) {
            const double ek = std::exp(_eigenVec[k] * t);
            for (int i = 0; i < N; ++i) {
                const double L_scaled = _L[i][k] * ek;
                // Direct pointer math for the inner loop to force contiguous vectorization
                double* row_tmp = &tmp[i * N];
                const double* row_R = _R[k]; 
                
                for (int j = 0; j < N; ++j) {
                    row_tmp[j] += L_scaled * row_R[j];
                }
            }
        }

        // Copy back to class member
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                _mat[i][j] = tmp[i * N + j];
            }
        }

        _cachedT = t;
    }

    replacementModel* _pb;

    double _eigenVec[N];
    double _L[N][N];
    double _R[N][N];

    mutable double  _mat[N][N];
    mutable MDOUBLE _cachedT;
};

#endif