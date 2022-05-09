#pragma once

class SymetricMatrix {

public:

    // Constructor
    SymetricMatrix(double c = 0)
    {
        for (int i = 0; i < 10; ++i)
        {
            m_m[i] = c;
        }
    }

    SymetricMatrix(double m11, double m12, double m13, double m14,
        double m22, double m23, double m24,
        double m33, double m34,
        double m44) 
    {
        m_m[0] = m11;  m_m[1] = m12;  m_m[2] = m13;  m_m[3] = m14;
        m_m[4] = m22;  m_m[5] = m23;  m_m[6] = m24;
        m_m[7] = m33;  m_m[8] = m34;
        m_m[9] = m44;
    }

    // Make plane
    SymetricMatrix(double a, double b, double c, double d)
    {
        m_m[0] = a * a;  m_m[1] = a * b;  m_m[2] = a * c;  m_m[3] = a * d;
        m_m[4] = b * b;  m_m[5] = b * c;  m_m[6] = b * d;
        m_m[7] = c * c; m_m[8] = c * d;
        m_m[9] = d * d;
    }

    double operator[](int c) const { return m_m[c]; }

    // Determinant
    double Det(int a11, int a12, int a13,
        int a21, int a22, int a23,
        int a31, int a32, int a33)
    {
        double det = m_m[a11] * m_m[a22] * m_m[a33] + m_m[a13] * m_m[a21] * m_m[a32] + m_m[a12] * m_m[a23] * m_m[a31]
            - m_m[a13] * m_m[a22] * m_m[a31] - m_m[a11] * m_m[a23] * m_m[a32] - m_m[a12] * m_m[a21] * m_m[a33];
        return det;
    }

    const SymetricMatrix operator+(const SymetricMatrix& n) const
    {
        return SymetricMatrix(m_m[0] + n[0], m_m[1] + n[1], m_m[2] + n[2], m_m[3] + n[3],
            m_m[4] + n[4], m_m[5] + n[5], m_m[6] + n[6],
            m_m[7] + n[7], m_m[8] + n[8],
            m_m[9] + n[9]);
    }

    SymetricMatrix& operator+=(const SymetricMatrix& n)
    {
        m_m[0] += n[0];   m_m[1] += n[1];   m_m[2] += n[2];   m_m[3] += n[3];
        m_m[4] += n[4];   m_m[5] += n[5];   m_m[6] += n[6];   m_m[7] += n[7];
        m_m[8] += n[8];   m_m[9] += n[9];
        return *this;
    }

    void Clear()
    {
        for (int i = 0; i < 10; ++i)
        {
            m_m[i] = 0;
        }
    }

private:
    double m_m[10];
};
