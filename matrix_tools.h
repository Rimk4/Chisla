#include <iostream>
#include <vector>

using namespace std;

#ifndef MATRIX_TOOLS
#define MATRIX_TOOLS

class ColumnVector
{
public:
	ColumnVector(const int n) : n(n), data(new double[n] {}) {}
	ColumnVector(const ColumnVector& v);

	ColumnVector(double* data, int n)
	{
		this->n = n;
		this->data = new double[this->n];
		for (size_t i = 0; i < this->n; i++)
		{
			this->data[i] = data[i];
		}
	}

	~ColumnVector();

	void setRand()
	{
		for (size_t i = 0; i < this->n; i++)
		{
			this->data[i] = (double)rand() / RAND_MAX;
		}
	}

	void print() const
	{
		for (size_t i = 0; i < this->n; i++)
		{
			cout << this->data[i] << " ";
		}
		cout << endl;
	}

	ColumnVector& operator=(const ColumnVector& v)
	{
		this->n = v.getN();
		for (size_t i = 0; i < this->n; i++)
		{
			this->data[i] = *(v[i]);
		}

		return *this;
	}

	ColumnVector operator+ (ColumnVector v)
	{
		ColumnVector result(v.getN());

		for (size_t i = 0; i < v.getN(); i++)
		{
			*(result[i]) = this->data[i] + *(v[i]);
		}

		return result;
	}

	double* operator[] (int index) const
	{
		return &(this->data[index]);
	}

	double operator* (const ColumnVector& v) const
	{
		double scalar = 0.0;

		for (size_t i = 0; i < this->n; i++)
		{
			scalar += this->data[i] * *(v[i]);
		}

		return scalar;
	}

	ColumnVector operator* (double alpha)
	{
		ColumnVector result(this->n);
		for (size_t i = 0; i < this->n; i++)
		{
			*(result[i]) = this->data[i] * alpha;
		}

		return result;
	}

	ColumnVector operator/ (double k) const
	{
		ColumnVector result(this->n);
		for (size_t i = 0; i < this->n; i++)
		{
			*(result[i]) = this->data[i] / k;
		}

		return result;
	}

	double getNorm() const
	{
		double norma = 0.0;

		for (size_t i = 0; i < this->n; i++)
		{
			norma += pow(this->data[i], 2);
		}

		return (double)sqrt(norma);
	}

	int getN() const
	{
		return this->n;
	}

private:
	double* data;
	int n;
};


class SquareMatrix
{
public:
	SquareMatrix(int n);
	SquareMatrix(double** data, int n);
	SquareMatrix(const SquareMatrix& matr_to_copy);
	SquareMatrix(vector<vector<double>> matr_to_copy)
	{
		this->n = matr_to_copy.size();
		this->data = new double* [this->n];
		for (size_t i = 0; i < this->n; i++)
		{
			this->data[i] = new double[this->n];
			for (size_t j = 0; j < this->n; j++)
			{
				this->data[i][j] = matr_to_copy[i][j];
			}
		}
	}
	//about constructors
	//https://cplusplus.com/articles/y8hv0pDG/#:~:text=the%20compiler%2Dprovided,heap%20corruption%20results.

	~SquareMatrix();

	void print() const;

	int getN() const
	{
		return this->n;
	}

	double** getData() const
	{
		return this->data;
	}

	void setRandValues()
	{
		for (size_t i = 0; i < this->n; i++)
		{
			for (size_t j = 0; j < this->n; j++)
			{
				this->data[i][j] = (double)rand() / RAND_MAX;
			}
		}
	}

	void identityMatrix()
	{
		for (size_t i = 0; i < this->n; i++)
		{
			this->data[i][i] = 1;
		}
	}

	void setRandDiag()
	{
		for (size_t i = 0; i < this->n; i++)
		{
			this->data[i][i] = (double)rand() / RAND_MAX;
		}
	}

	SquareMatrix getDiag()
	{
		SquareMatrix D(n);

		for (size_t i = 0; i < n; i++)
		{
			D[i][i] = data[i][i];
		}

		return D;
	}

	SquareMatrix T();


	SquareMatrix Inverse()
	{
		SquareMatrix A(this->data, this->n);
		double temp;


		SquareMatrix E(this->n);
		E.identityMatrix();

		for (int k = 0; k < this->n; k++)
		{
			temp = A[k][k];

			for (int j = 0; j < this->n; j++)
			{
				A[k][j] /= temp;
				E[k][j] /= temp;
			}

			for (int i = k + 1; i < this->n; i++)
			{
				temp = A[i][k];

				for (int j = 0; j < this->n; j++)
				{
					A[i][j] -= A[k][j] * temp;
					E[i][j] -= E[k][j] * temp;
				}
			}
		}

		for (int k = this->n - 1; k > 0; k--)
		{
			for (int i = k - 1; i >= 0; i--)
			{
				temp = A[i][k];

				for (int j = 0; j < this->n; j++)
				{
					A[i][j] -= A[k][j] * temp;
					E[i][j] -= E[k][j] * temp;
				}
			}
		}

		for (int i = 0; i < this->n; i++)
			for (int j = 0; j < this->n; j++)
				A[i][j] = E[i][j];

		return A;
	}

	double* operator[] (int index) const
	{
		return data[index];
	}

	SquareMatrix operator+ (SquareMatrix matr_to_sum);
	SquareMatrix operator* (SquareMatrix matr);
	SquareMatrix operator* (const double k);
	vector<double> operator* (vector<double> v);

	SquareMatrix& operator=(const SquareMatrix& m); //https://cplusplus.com/articles/y8hv0pDG/#:~:text=When%20do%20I%20need%20to%20write%20an,need%20to%20write%20a%20custom%20assignment%20operator.

	ColumnVector operator* (ColumnVector& vec)
	{
		ColumnVector result(vec.getN());

		for (size_t i = 0; i < this->n; i++)
		{
			for (size_t j = 0; j < this->n; j++)
			{
				*(result[i]) += *(vec[j]) * this->data[i][j];
			}
		}


		return result;
	}

private:
	int n;
	double** data;
};

#endif