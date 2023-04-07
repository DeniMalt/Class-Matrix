#pragma once
#include<iostream>
#include<fstream>
#include<ostream>
using namespace std;

template <typename T>
class Matrix {
private:
    T* pMatrix;
    int m;
    int n;
public:
    Matrix(int Rows, int Cols) 
    {
        m = Rows;
        n = Cols;
        pMatrix = new T[m * n];
    }

    Matrix(int Size) 
    {
        m = Size;
        n = Size;
        pMatrix = new T[m * n];
    }

    Matrix(const Matrix& matrix) 
    {
        m = matrix.m;
        n = matrix.n;
        pMatrix = new T[m * n];

        for (int i = 0; i < m * n; i++) 
        {
            pMatrix[i] = matrix.pMatrix[i];
        }
    }

    ~Matrix() 
    {
        delete[] pMatrix;
    }

    Matrix& operator = (const Matrix& matrix) 
    {
        this->~Matrix();
        m = matrix.m;
        n = matrix.n;
        pMatrix = new T[n * m];
        for (int i = 0; i < n * m; i++) 
        {
            pMatrix[i] = matrix.pMatrix[i];
        }
        return *this;
    }

    bool operator == (const Matrix& matrix) 
    {
        if (m == matrix.m && n == matrix.n) 
        {
            for (int i = 0; i < m * n; i++) 
            {
                if (pMatrix[i] != matrix.pMatrix[i]) 
                {
                    return false;
                }
            }
            return true;
        }
    }

    bool operator != (const Matrix& matrix) 
    {
        if (m == matrix.m && n == matrix.n) 
        {
            for (int i = 0; i < m * n; i++) 
            {
                if (pMatrix[i] != matrix.pMatrix[i]) 
                {
                    return true;
                }
            }
            return false;
        }
    }

    Matrix operator+ (const Matrix& matrix) 
    {
        if (m == matrix.m && n == matrix.n) 
        {
            Matrix result(m, n);
            for (int i = 0; i < m * n; i++) 
            {
                result.pMatrix[i] = pMatrix[i] + matrix.pMatrix[i];
            }
            return result;
        }
    }

    Matrix operator- (const Matrix& matrix)
    {
        if (m == matrix.m && n == matrix.n) 
        {
            Matrix result(m, n);
            for (int i = 0; i < m * n; i++) 
            {
                result.pMatrix[i] = pMatrix[i] - matrix.pMatrix[i];
            }
            return result;
        }
    }

    Matrix operator* (const Matrix& matrix) 
    {
        if (n == matrix.m)
        {
            Matrix result(m, matrix.n);
            for (int i = 0; i < m; i++) 
            {
                for (int j = 0; j < matrix.n; j++) 
                {
                    result.pMatrix[i * result.n + j] = 0;
                    for (int k = 0; k < n; k++) 
                    {
                        result.pMatrix[i * result.n + j] += pMatrix[i * n + k] * matrix.pMatrix[k * matrix.n + j];
                    }
                }
            }
            return result;
        }
    }


    Matrix operator* (const int k)
    {
        Matrix result(n, m);
        for (int i = 0; i < m; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                result.pMatrix[j * m + i] = pMatrix[j * n + i] * k;
            }
        }
        return result;
    }


    Matrix operator/ (const int k)
    {
        Matrix result(n, m);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result.pMatrix[j * m + i] = pMatrix[j * n + i] / k;
            }
        }
        return result;
    }

    Matrix operator* (const double k)
    {
        Matrix result(n, m);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result.pMatrix[j * m + i] = pMatrix[j * n + i] * k;
            }
        }
        return result;
    }

    Matrix operator/ (const double k)
    {
        Matrix result(n, m);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result.pMatrix[j * m + i] = pMatrix[j * n + i] / k;
            }
        }
        return result;
    }
    
    void operator*= (const Matrix& matrix)
    {
        *this = *this * matrix;
    }

    Matrix pow(const int k)
    {
        if (m == n)
        {
            if (k == 0)
            {
                return this->Matrix_E();
            }
            if (k > 0)
            {
                Matrix result = *this;
                for (int i = 1; i < k; i++)
                {
                    result *= *this;
                }
                return result;
            }
            if (k < 0)
            {
                Matrix result = InverseMatrix();
                for (int i = 1; i < abs(k); i++)
                {
                    result *= InverseMatrix();
                }
                return result;

            }
        }
        else
        {
            cout << "Error" << endl;
            return 0;
        }
    }

    Matrix Matrix_E()
    {
        Matrix result(m, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (i == j)
                {
                    result.pMatrix[j * n + i] = 1;
                }
                else
                {
                    result.pMatrix[j * n + i] = 0;
                }
            }
        }
        return result;
    }


    Matrix Matrix_of_minors(int rows, int cols)
    {
        if (rows <= m && cols <= n)
        {
            Matrix result(m - 1, n - 1);
            int curi = 0;
            int curj = 0;
            for (int i = 0; i < m; i++)
            {
                if (i == rows - 1)
                    curi--;
                for (int j = 0; j < n; j++)
                {
                    if (j == cols - 1)
                        curj--;
                    if (i != rows - 1 && j != cols - 1)
                        result.pMatrix[curi * (n - 1) + curj] = pMatrix[i * n + j];
                    curj++;
                }
                curj = 0;
                curi++;
            }
            return result;
        }
    }

 
    Matrix Transponse()
    {
        Matrix result(n, m);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result.pMatrix[j * m + i] = pMatrix[i * n + j];
            }
        }
        return result;
    }
   
    T Minor(int Row, int Col)
    {
        if (m == n && Row <= m && Col <= n)
        {
            return Matrix_of_minors(Row, Col).Det();
        }
    }
    
    T AlgCompl(int Row, int Col)
    {
        if (m == n && Row <= m && Col <= n)
        {
            return (((Row + Col) % 2 == 0) ? 1 : -1) * Minor(Row, Col);
        }
    }

    T Det()
    {
        if (m == n)
        {
            if (m == 1)
            {
                return pMatrix[0];
            }

            else if (m == 2)
            {
                return pMatrix[0] * pMatrix[3] - pMatrix[1] * pMatrix[2];
            }

            else
            {
                T det = pMatrix[0] * AlgCompl(1, 1);
                for (int i = 2; i <= m; i++)
                {
                    det += pMatrix[(i - 1) * n] * AlgCompl(i, 1);
                }
                return det;
            }
        }
    }
    
    Matrix AddNumberStr(int a, int Row)
    {
        Matrix result(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (Row == j)
                {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i] + a;
                }

                else {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i];
                }
            }
        }
        return result;
    }

    Matrix AddNumberStr(double a, int Row)
    {
        Matrix result(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (Row == j)
                {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i] + a;
                }

                else {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i];
                }
            }
        }
        return result;
    }

    Matrix MultiNumberStr(int a, int Row)
    {
        Matrix result(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (Row == j)
                {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i] * a;
                }

                else {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i];
                }
            }
        }
        return result;
    }

    Matrix MultiNumberStr(double a, int Row)
    {
        Matrix result(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (Row == j)
                {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i] * a;
                }

                else {
                    result.pMatrix[j * m + i] = pMatrix[j * n + i];
                }
            }
        }
        return result;
    }

    Matrix Diffstrs(int st1, int st2)
    {
        int str1 = st1 - 1;
        int str2 = st2 - 1;
        if (str1 != str2)
        {
            Matrix result(m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (str1 = j)
                    {
                        result.pMatrix[j * m + i] = pMatrix[j * m + i] - pMatrix[str2 * n + i];
                        result.pMatrix[str2 * m + i] = pMatrix[str2 * m + i];
                    }
                    else
                    {
                        result.pMatrix[j * m + i] = pMatrix[j * m + i];
                    }
                }
            }
            return result;
        }
        else
        {
            cout << "Error" << endl;
            return 0;
        }
    }

    Matrix Summstrs(int st1, int st2)
    {
        int str1 = st1 - 1;
        int str2 = st2 - 1;
        if (str1 != str2)
        {
            Matrix result(m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (str1 == j)
                    {
                        result.pMatrix[j * m + i] = pMatrix[j * m + i] + pMatrix[str2 * n + i];
                        result.pMatrix[str2 * m + i] = pMatrix[str2 * m + i];
                    }
                    else
                    {
                        result.pMatrix[j * m + i] = pMatrix[j * m + i];
                    }
                }
            }
            return result;
        }
        else
        {
            cout << "Error" << endl;
            return 0;
        }
    }

    Matrix InverseMatrix()
    {   

        if (m == n && Det() != 0)
        {
            if (m == 1 && n == 1)
            {
                Matrix result(m, n);
                result.pMatrix[0] = 1 / pMatrix[0];
                return result;
            }
            else
            {
                Matrix result(m, n);
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        result.pMatrix[i * m + j] = this->AlgCompl(j + 1, i + 1) / Det();
                    }
                }
                return result;
            }
        }
        else
        {
            cout << "Error" << endl;
        }
    }

    Matrix ReplaceStrs(int st1, int st2)
    {
        int str1 = st1 - 1;
        int str2 = st2 - 1;
        if (str1 != str2)
        {
            Matrix result(m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {   
                    if (str2 == j)
                    {
                        result.pMatrix[str2 * m + i] = pMatrix[str1 * m + i];
                        result.pMatrix[str1 * m + i] = pMatrix[str2 * m + i];
                    }
                    else
                    {
                        result.pMatrix[j * m + i] = pMatrix[j * m + i];
                    }
                }
            }
            return result;
        }
        else
        {
            cout << "Error" << endl;
            return 0;
        }
    }

    Matrix ReplaceCol(int c1, int c2)
    {
        int col1 = c1 - 1;
        int col2 = c2 - 1;
        if (col1 != col2)
        {
            Matrix result(m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (col2 == i)
                    {
                        result.pMatrix[j * m + col1] = pMatrix[j * m + col2];
                        result.pMatrix[j * m + col2] = pMatrix[j * m + col1];
                    }
                    else
                    {
                        result.pMatrix[j * m + i] = pMatrix[j * m + i];
                    }
                }
            }
            return result;
        }
        else
        {
            cout << "Error" << endl;
            return 0;
        }
    }

    friend istream& operator>> <> (istream& s, Matrix& matrix);
    friend ostream& operator<< <> (ostream& s, const Matrix& matrix);
};


template <typename T>
istream& operator >> (istream& s, Matrix<T>& matrix) 
{
    for (int i = 0; i < matrix.m * matrix.n; i++) 
    {
        s >> matrix.pMatrix[i];
    }
    return s;
}

template <typename T>
ostream& operator << (ostream& s, const Matrix<T>& matrix) 
{
    for (int i = 0; i < matrix.m; i++) 
    {
        for (int j = 0; j < matrix.n; j++) 
        {
            s << matrix.pMatrix[i * matrix.n + j] << " ";
        }
        s << endl;
    }
    return s;
}




