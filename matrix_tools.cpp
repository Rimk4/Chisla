#include <iostream>
#include <vector>
#include "matrix_tools.h"

using namespace std;

ColumnVector::ColumnVector(const ColumnVector& v) : n(v.getN())
{
	this->data = new double[n] {};
	for (size_t i = 0; i < this->n; i++)
	{
		this->data[i] = *(v[i]);
	}
}

ColumnVector::~ColumnVector()
{
	delete[] data;
}

void SquareMatrix::print() const
{
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			cout << this->data[i][j] << " ";
		}
		cout << endl;
	}
}

SquareMatrix SquareMatrix::T()
{
	SquareMatrix result(this->n);

	for (size_t i = 0; i < this->n; i++)
	{
		for (size_t j = 0; j < this->n; j++)
		{
			result[i][j] = this->data[j][i];
		}
	}

	return result;
}


//конструктортдля генерации "нулевой" матрицы
SquareMatrix::SquareMatrix(int n)
{
	this->n = n;
	this->data = new double* [n];
	for (size_t i = 0; i < n; i++)
	{
		this->data[i] = new double[n] {};
	}
}

//конструктор для копирования от массива
SquareMatrix::SquareMatrix(double** data, int n)
{
	this->n = n;

	this->data = new double* [this->n];
	for (size_t i = 0; i < this->n; i++)
	{
		this->data[i] = new double[this->n] {};
		for (size_t j = 0; j < this->n; j++)
		{
			this->data[i][j] = data[i][j];
		}
	}
}

//Конструктор копирования - deep copy (copy what is pointed to).
SquareMatrix::SquareMatrix(const SquareMatrix& matr_to_copy)
{
	this->n = matr_to_copy.getN();

	this->data = new double* [this->n];
	for (size_t i = 0; i < this->n; i++)
	{
		this->data[i] = new double[this->n] {};
		for (size_t j = 0; j < this->n; j++)
		{
			this->data[i][j] = matr_to_copy[i][j];
		}
	}
}

//деструктор
SquareMatrix::~SquareMatrix()
{
	for (size_t i = 0; i < this->n; i++)
	{
		delete[] this->data[i];
	}
	delete[] this->data;
}

//оператор сложения матриц
SquareMatrix SquareMatrix::operator+ (SquareMatrix matr_to_sum)
{
	SquareMatrix result(matr_to_sum.getN());

	for (size_t i = 0; i < result.getN(); i++)
	{
		for (size_t j = 0; j < result.getN(); j++)
		{
			result[i][j] = this->data[i][j] + matr_to_sum[i][j];
		}
	}

	return result;
}

//оператор умножения матриц
SquareMatrix SquareMatrix::operator* (SquareMatrix matr)
{
	double** result = new double* [this->n];
	for (size_t i = 0; i < this->n; i++)
	{
		result[i] = new double[this->n] {};
	}

	for (size_t i = 0; i < this->n; i++)
	{
		for (size_t j = 0; j < this->n; j++)
		{
			for (size_t k = 0; k < this->n; k++)
			{
				result[i][j] += this->data[i][k] * matr[k][j];
			}
		}
	}

	return SquareMatrix(result, this->n);
}

//оператор присваивания
SquareMatrix& SquareMatrix::operator=(const SquareMatrix& m)
{
	this->n = m.getN();

	for (size_t i = 0; i < this->n; i++)
	{
		for (size_t j = 0; j < this->n; j++)
		{
			this->data[i][j] = m[i][j];
		}
	}

	return *this;
}

//оператор умножения на число
SquareMatrix SquareMatrix::operator* (const double k)
{
	SquareMatrix result(this->n);

	for (size_t i = 0; i < this->n; i++)
	{
		for (size_t j = 0; j < this->n; j++)
		{
			result[i][j] = this->data[i][j] * k;
		}
	}

	return result;
}

//оператор умножения на столбец
vector<double> SquareMatrix::operator* (vector<double> v)
{
	vector<double> result(this->n, 0);

	for (size_t i = 0; i < this->n; i++)
	{
		for (size_t j = 0; j < this->n; j++)
		{
			result[i] += this->data[i][j] * v[j];
		}
	}

	return result;
}