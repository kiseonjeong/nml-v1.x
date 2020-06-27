#include "stdafx.h"
#include "Decomposition.h"

namespace nml
{
	const numat decomposition::operator[](const int idx) const
	{
		// Check the index
		assert(idx >= 0 && idx < 2);

		return D[idx];
	}

	void decomposition::checkEmpty(const numat& mat) const
	{
		// Check the matrix is empty or not
		assert(mat.empty() == false);
	}

	void decomposition::checkSquareness(const numat& mat) const
	{
		// Check the matrix is square or not
		assert(mat.rows == mat.cols);
	}

	const vector<numat> decomposition::result() const
	{
		// Get the decomposition result
		return D;
	}

	LUD::LUD()
	{
		// Set an object
		setObject();
	}

	LUD::LUD(const numat& mat)
	{
		// Set an object
		setObject();

		// Set a partial pivot as a default mode
		pivot(LU_PARTIAL_PIVOT);

		// Decompose the matrix
		decomp(mat);
	}

	LUD::LUD(const numat& mat, const pivotMode mode)
	{
		// Set an object
		setObject();

		// Set a pivot mode
		pivot(mode);

		// Decompose the matrix
		decomp(mat);
	}

	LUD::LUD(const LUD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	LUD::~LUD()
	{
		// Clear the object
		clearObject();
	}

	LUD& LUD::operator=(const LUD& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void LUD::setObject()
	{
		// Set the parameters
		setType(*this);
		mode = LU_PARTIAL_PIVOT;

		// Set the memories
		P.clear();
		M.clear();
		D.clear();
	}

	void LUD::copyObject(const object& obj)
	{
		// Do down casting
		LUD* _obj = (LUD*)&obj;

		// Copy the parameters
		mode = _obj->mode;

		// Copy the memories
		P = _obj->P;
		M = _obj->M;
		D = _obj->D;
	}

	void LUD::clearObject()
	{
		// Clear the memories
		P.clear();
		M.clear();
		D.clear();
	}

	void LUD::pivot(const pivotMode mode)
	{
		// Set a pivot mode
		this->mode = mode;
	}

	const vector<numat> LUD::decomp(const numat& mat)
	{
		// Check the matrix
		checkEmpty(mat);

		// Check the matrix
		checkSquareness(mat);

		// Clear the old memories
		P.clear();
		M.clear();
		D.clear();

		// Copy the input matrix
		numat A = mat;

		// Factorize the matrix
		factorize(mat);

		// Create a triangle matrix
		vector<numat> triangle;
		for (int i = mat.rows - 1; i >= 0; i--)
		{
			triangle.push_back(M[i] * P[i]);
		}

		// Create a decomposition matrix vector
		numat temp = triangle[0];
		for (int i = 1; i < (int)triangle.size(); i++)
		{
			temp *= triangle[i];
		}
		D.push_back(temp.I());
		D.push_back(temp * A);

		return D;
	}

	void LUD::factorize(const numat& mat)
	{
		// Factorize the matrix with pivoting
		numat temp = mat;
		for (int i = 0; i < temp.cols; i++)
		{
			// Find the row and column maximum index
			int rowMaxIdx = i;
			int colMaxIdx = i;
			double maxVal = temp(rowMaxIdx, colMaxIdx);
			for (int j = i; j < temp.rows; j++)
			{
				// Check the pivot mode
				if (mode == LU_PARTIAL_PIVOT)
				{
					if (temp(j, i) > maxVal)
					{
						maxVal = temp(j, i);
						rowMaxIdx = j;
						colMaxIdx = i;
					}
				}
				else
				{
					for (int k = i; k < temp.cols; k++)
					{
						if (temp(j, k) > maxVal)
						{
							maxVal = temp(j, k);
							rowMaxIdx = j;
							colMaxIdx = k;
						}
					}
				}
			}

			// Swap the matrix
			temp.swapr(i, rowMaxIdx);
			temp.swapc(i, colMaxIdx);

			// Save the result
			P.push_back(numat::permut(temp.rows, i, rowMaxIdx));
			M.push_back(numat::eyes(temp.rows));
			for (int j = i + 1; j < temp.rows; j++)
			{
				M[M.size() - 1](j, i) = -temp(j, i) / temp(i, i);
				for (int k = 0; k < temp.cols; k++)
				{
					temp(j, k) += temp(i, k) * M[M.size() - 1](j, i);
				}
			}
		}
	}

	QRD::QRD()
	{
		// Set an object
		setObject();
	}

	QRD::QRD(const numat& mat)
	{
		// Set an object
		setObject();

		// Decompose the matrix
		decomp(mat);
	}

	QRD::QRD(const QRD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	QRD::~QRD()
	{
		// Clear the object
		clearObject();
	}

	QRD& QRD::operator=(const QRD& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void QRD::setObject()
	{
		// Set the parameters
		setType(*this);

		// Set the memories
		T.clear();
		D.clear();
	}

	void QRD::copyObject(const object& obj)
	{
		// Do down casting
		QRD* _obj = (QRD*)&obj;

		// Copy the memories
		T = _obj->T;
		D = _obj->D;
	}

	void QRD::clearObject()
	{
		// Clear the memories
		T.clear();
		D.clear();
	}

	const vector<numat> QRD::decomp(const numat& mat)
	{
		// Check Matrix
		checkEmpty(mat);

		// Clear the old memories
		T.clear();
		D.clear();

		// Check the row information of the matrix
		if (mat.rows == 1)
		{
			D.push_back(numat::ones(msize(1)));
			D.push_back(mat);

			return D;
		}

		// Copy the input matrix
		numat temp = mat;
		numat A = mat;

		// Do House Holder Transform on the matrix
		for (int i = 0; i < A.rows - 1; i++)
		{
			int j = (i != 0) ? 1 : 0;
			if (j - temp.cols == 0)
			{
				break;
			}
			numat S = reduction(temp, j);
			numat Q = houseHolder(S);
			temp = Q * S;
			T.push_back(extension(Q, i, A.rows));
		}

		// Create a decomposition matrix vector
		const int count = (int)T.size();
		D.push_back(T[0].T());
		for (int i = 1; i < count; i++)
		{
			D[0] *= T[i].T();
		}
		D.push_back(T[count - 1]);
		for (int i = count - 2; i >= 0; i--)
		{
			D[1] *= T[i];
		}
		D[1] *= A;

		return D;
	}

	const numat QRD::reduction(const numat& mat, int flag) const
	{
		// Check the reduction flag
		if (flag == 0)
		{
			// Get the input matrix
			return mat;
		}
		else
		{
			// Get a reduction matrix
			return mat.minor(0, 0);
		}
	}

	const numat QRD::houseHolder(const numat& mat) const
	{
		// Do House Holder Transform on the matrix
		numat a = mat.col(0);
		numat e = reflection(a, mat.rows);
		numat u = a - e;
		numat v = u / 2.0;
		numat vvt = v * v.T();
		
		return numat::eyes(mat.rows) - (vvt * (2.0 / vecsize(vvt)));
	}

	const numat QRD::reflection(const numat& mat, int length) const
	{
		// Create a reflection matrix
		numat result(msize(length, 1), 0.0);
		result(0) = vecsize(mat);

		return result;
	}

	const double QRD::vecsize(const numat& mat) const
	{
		// Calculate a size of the input vector
		return std::sqrt(numat::sum(mat.mul(mat)));
	}

	const numat QRD::extension(const numat& mat, const int offset, const int length) const
	{
		// Get an extension matrix
		numat result(msize(length, length), 0.0);
		for (int i = 0; i < offset; i++)
		{
			result(i, i) = 1.0;
		}
		for (int i = 0; i < mat.rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				result(offset + i, offset + j) = mat(i, j);
			}
		}

		return result;
	}

	EVD::EVD()
	{
		// Set an object
		setObject();
	}

	EVD::EVD(const numat& mat)
	{
		// Set an object
		setObject();

		// Decompose the matrix
		decomp(mat);
	}

	EVD::EVD(const EVD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	EVD::~EVD()
	{
		// Clear the object
		clearObject();
	}

	EVD& EVD::operator=(const EVD& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void EVD::setObject()
	{
		// Set the parameters
		setType(*this);

		// Set the memories
		D.clear();
	}

	void EVD::copyObject(const object& obj)
	{
		// Do down casting
		EVD* _obj = (EVD*)&obj;

		// Copy the memories
		D = _obj->D;
	}

	void EVD::clearObject()
	{
		// Clear the memories
		D.clear();
	}

	const vector<numat> EVD::decomp(const numat& mat)
	{
		// Check the matrix
		checkEmpty(mat);

		// Check the matrix
		checkSquareness(mat);

		// Clear the old memories
		D.clear();

		// Initialize the parameters
		const int maxIter = 100 * max(mat.rows, mat.cols);
		int loopIter = 0;
		bool convFlag = false;

		// Do eigen value decomposition using QRD
		QRD qr(mat);
		numat eigVec = qr[0];
		numat eigVal = qr[1] * qr[0];
		numat tempVal = eigVal;
		while (true)
		{
			// Do QR decomposition
			qr.decomp(eigVal);

			// Calculate an eigen value and an eigen vector
			eigVal = qr[1] * qr[0];
			convFlag = isConvergence(eigVal, tempVal);
			eigVec *= qr[0];

			// Check the convergence flag
			if (convFlag == true)
			{
				break;
			}
			else
			{
				// Check the stop condition
				if (loopIter == maxIter)
				{
					break;
				}
				else
				{
					tempVal = eigVal;
					loopIter++;
				}
			}
		}

		// Create a decomposition matrix vector
		D.push_back(eigVec);
		D.push_back(eigVal);

		return D;
	}

	const bool EVD::isConvergence(const numat& curr, const numat& prev) const
	{
		// Calculate a delta between matrices
		numat delta = numat::sqrt(curr.mul(curr) - prev.mul(prev));

		// Check the matrix is convergence or not
		for (int i = 0; i < delta.length(); i++)
		{
			if (numat::almost(delta(i), 0.0, 5) == false)
			{
				return false;
			}
		}

		return true;
	}

	SVD::SVD()
	{
		// Set an object
		setObject();
	}

	SVD::SVD(const numat& mat)
	{
		// Set an object
		setObject();

		// Decompose the matrix
		decomp(mat);
	}

	SVD::SVD(const SVD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	SVD::~SVD()
	{
		// Clear the object
		clearObject();
	}

	SVD& SVD::operator=(const SVD& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	const numat SVD::operator[](const int idx) const
	{
		// Check the index
		assert(idx >= 0 && idx < 3);

		return D[idx];
	}

	void SVD::setObject()
	{
		// Set the parameters
		setType(*this);
	}

	const vector<numat> SVD::decomp(const numat& mat)
	{
		// Check the matrix
		checkEmpty(mat);

		// Clear the old memories
		D.clear();

		// Initialize the parameters
		const int maxIter = 100 * max(mat.rows, mat.cols);
		int loopIter = 0;

		// Do eigen value decomposition using QRD
		QRD qru(mat);
		QRD qrv(qru[1].T());
		numat T = qrv[1];
		numat U = qru[0];
		numat S = qrv[1];
		numat V = qrv[0];
		while (true)
		{
			// Do QR decomposition
			qru.decomp(S.T());
			qrv.decomp(qru[1].T());
			U *= qru[0];
			V *= qrv[0];
			S = qrv[1];

			// Check the convergence flag
			if (isConvergence(S, T) == true)
			{
				break;
			}
			else
			{
				// Check the stop condition
				if (loopIter == maxIter)
				{
					break;
				}
				else
				{
					T = S;
					loopIter++;
				}
			}
		}

		// Create a decomposition matrix vector
		S = S.T();
		int k = (mat.cols > mat.rows) ? mat.rows : mat.cols;
		int l = mat.rows;
		for (int i = 0; i < k; i++)
		{
			if (S(i, i) < 0.0)
			{
				S(i, i) = -S(i, i);
				for (int j = 0; j < l; j++)
				{
					U(j, i) = -U(j, i);
				}
			}
		}
		D.push_back(U);
		D.push_back(S);
		D.push_back(V);

		return D;
	}

	numat SVD::pinv(const numat& mat)
	{
		// Factorize the matrix
		SVD sv(mat);

		// Inverse the elements on the eigen matrix
		numat E = sv[1];
		for (int i = 0; i < E.rows; i++)
		{
			for (int j = 0; j < E.cols; j++)
			{
				if (numat::almost(E(i, j), 0.0, 5) == false)
				{
					E(i, j) = 1.0 / E(i, j);
				}
			}
		}

		return sv[2] * E.T() * sv[0].T();
	}
}