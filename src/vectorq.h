#pragma once

#include <Eigen/Core>
#include <gmpxx.h>

template <>
struct Eigen::NumTraits<mpq_class> : Eigen::GenericNumTraits<mpq_class> {
    typedef mpq_class Real;
    typedef mpq_class NonInteger;
    typedef mpq_class Nested;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }
    static inline int digits10() { return 0; }

    enum {
        IsInteger = 0,
        IsSigned = 1,
        IsComplex = 0,
        RequireInitialization = 1,
        ReadCost = 6,
        AddCost = 150,
        MulCost = 100
    };
};

typedef Eigen::Matrix<mpq_class, 3, 1> Vector3q;
typedef Eigen::Matrix<mpq_class, 3, 3> Matrix3q;
