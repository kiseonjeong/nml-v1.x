#include "stdafx.h"
#include "Numat.h"

namespace nml
{
	numat::numat()
	{
		// Set an object
		setObject();
	}

	numat::numat(const msize& size)
	{
		// Set an object
		setObject();

		// Create a matrix
		create(size);
	}

	numat::numat(const msize& size, double val)
	{
		// Set an object
		setObject();

		// Create a matrix
		create(size, val);
	}

	numat::numat(const numem<double>& data)
	{
		// Set an object
		setObject();

		// Create a matrix
		copy(data);
	}

	numat::numat(const numat& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	numat::~numat()
	{
		// Clear the object
		clearObject();
	}

	numat& numat::operator=(const numat& mat)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(mat);

		return *this;
	}

	void numat::setObject()
	{
		// Set the parameters
		setType(*this);
	}

	const numat numat::operator+(const numat& mat) const
	{
		return add(mat);
	}

	const numat numat::operator+(double val) const
	{
		return add(val);
	}

	numat& numat::operator+=(const numat& mat)
	{
		*this = add(mat);

		return *this;
	}

	numat& numat::operator+=(double val)
	{
		*this = add(val);

		return *this;
	}

	const numat numat::operator-(const numat& mat) const
	{
		return sub(mat);
	}

	const numat numat::operator-(double val) const
	{
		return sub(val);
	}

	numat& numat::operator-=(const numat& mat)
	{
		*this = sub(mat);

		return *this;
	}

	numat& numat::operator-=(double val)
	{
		*this = sub(val);

		return *this;
	}

	const numat numat::operator*(const numat& mat) const
	{
		return dot(mat);
	}

	const numat numat::operator*(double val) const
	{
		return mul(val);
	}

	numat& numat::operator*=(const numat& mat)
	{
		*this = dot(mat);

		return *this;
	}

	numat& numat::operator*=(double val)
	{
		*this = mul(val);

		return *this;
	}

	const numat numat::operator/(const numat& mat) const
	{
		return div(mat);
	}

	const numat numat::operator/(double val) const
	{
		return div(val);
	}

	numat& numat::operator/=(const numat& mat)
	{
		*this = div(mat);

		return *this;
	}

	numat& numat::operator/=(double val)
	{
		*this = div(val);

		return *this;
	}

	void numat::cout() const
	{
		// Check the data
		if (empty() == true)
		{
			std::cout << "[]" << endl;
		}
		else
		{
			if (rows == 1)
			{
				for (int i = 0; i < cols; i++)
				{
					if (i == 0)
					{
						std::cout << "[";
					}
					std::cout << (*this)(i);
					if (i == cols - 1)
					{
						std::cout << "]";
					}
					else
					{
						std::cout << ", ";
					}
				}
			}
			else if (cols == 1)
			{
				for (int i = 0; i < rows; i++)
				{
					if (i == 0)
					{
						std::cout << "[";
					}
					std::cout << (*this)(i);
					if (i == rows - 1)
					{
						std::cout << "]";
					}
					else
					{
						std::cout << ", " << endl;
					}
				}
			}
			else
			{
				for (int i = 0; i < rows; i++)
				{
					if (i == 0)
					{
						std::cout << "[";
					}
					for (int j = 0; j < cols; j++)
					{
						if (j == 0)
						{
							std::cout << "[";
						}
						std::cout << (*this)(i * cols + j);
						if (j == cols - 1)
						{
							std::cout << "]";
						}
						else
						{
							std::cout << ", ";
						}
					}
					if (i == rows - 1)
					{
						std::cout << "]";
					}
					else
					{
						std::cout << ", " << endl << " ";
					}
				}
			}
			std::cout << endl;
		}
		std::cout << endl;
	}

	const numem<double> numat::to_data() const
	{
		// Check the matrix
		if (empty() == true)
		{
			return numem<double>();
		}

		// Convert the matrix to the 'numem' data
		return numem<double>(*this);
	}

	numat numat::to_mat(const numem<double>& data)
	{
		// Check the data
		if (data.empty() == true)
		{
			return numat();
		}

		// Convert the data to the matrix
		return numat(data);
	}

	numat numat::to_mat(const vector<double>& vec)
	{
		// Check the vector
		if (vec.empty() == true)
		{
			return numat();
		}

		// Create a matrix
		numat result(msize((int)vec.size(), 1));

		// Initialize the matrix
		for (int i = 0; i < result.length(); i++)
		{
			result(i) = vec[i];
		}

		return result;
	}

	numat numat::to_mat(const vector<double>& vec, int nrows)
	{
		// Check the vector
		if (vec.empty() == true)
		{
			return numat();
		}

		// Convert the vector to the matrix
		numat temp = to_mat(vec);

		// Reshape the matrix
		return temp.reshape(nrows);
	}

	const numat numat::clone() const
	{
		// Check the matrix
		if (empty() == true)
		{
			return numat();
		}

		// Clone the matrix
		return numat(*this);
	}

	const numat numat::reshape(int nrows) const
	{
		// Check the matrix
		if (empty() == true)
		{
			return numat();
		}

		// Check the size
		assert(nrows > 0 && nrows <= length());

		// Check the dimension
		assert(length() % nrows == 0);

		// Check the current row information
		if (rows == nrows)
		{
			return numat(*this);
		}

		// Set the rows and the columns
		int rows = nrows;
		int cols = (int)ceil((double)length() / nrows);

		// Create a reshaped matrix
		numat result(msize(rows, cols), 0.0);
		for (int i = 0; i < length(); i++)
		{
			result(i) = (*this)(i);
		}

		return result;
	}

	numat numat::vappend(const vector<numat>& arr)
	{
		// Check the array
		assert(arr.empty() == false);
		const int mats = (int)arr.size();

		// Check the size
		vector<int> rinfo;
		int rows = 0;
		int cols = arr[0].cols;
		for (int i = 0; i < mats; i++)
		{
			if (i != 0)
			{
				assert(arr[0].cols == arr[i].cols);
			}
			rinfo.push_back(arr[i].rows);
			rows += arr[i].rows;
		}

		// Append the matrices
		numat result(msize(rows, cols));
		for (int i = 0, j = 0; i < mats; i++)
		{
			for (int k = 0; k < rinfo[i]; k++)
			{
				for (int l = 0; l < cols; l++)
				{
					result(j + k, l) = arr[i](k, l);
				}
			}
			j += rinfo[i];
		}

		return result;
	}

	const numat numat::vappend(const numat& mat) const
	{
		// Create an array
		vector<numat> arr;
		arr.push_back(*this);
		arr.push_back(mat);

		// Append the matrices
		return vappend(arr);
	}

	numat numat::happend(const vector<numat>& arr)
	{
		// Check the array
		assert(arr.empty() == false);
		const int mats = (int)arr.size();

		// Check the size
		vector<int> cinfo;
		int rows = arr[0].rows;
		int cols = 0;
		for (int i = 0; i < mats; i++)
		{
			if (i != 0)
			{
				assert(arr[0].rows == arr[i].rows);
			}
			cinfo.push_back(arr[i].cols);
			cols += arr[i].cols;
		}

		// Append the matrices
		numat result(msize(rows, cols));
		for (int i = 0, j = 0; i < mats; i++)
		{
			for (int k = 0; k < cinfo[i]; k++)
			{
				for (int l = 0; l < rows; l++)
				{
					result(l, j + k) = arr[i](l, k);
				}
			}
			j += cinfo[i];
		}

		return result;
	}

	const numat numat::happend(const numat& mat) const
	{
		// Create an array
		vector<numat> arr;
		arr.push_back(*this);
		arr.push_back(mat);

		// Append the matrices
		return happend(arr);
	}

	const numat numat::row(int idx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(idx >= 0 && idx < rows);

		// Get the selected row matrix
		numat result(msize(1, cols));
		for (int i = 0; i < cols; i++)
		{
			result(i) = (*this)(idx, i);
		}

		return result;
	}

	const numat numat::row(int bi, int ei, int ii) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(bi >= 0 && bi < rows);
		assert(ei >= 0 && ei < rows);
		assert(bi < ei);
		if (ii < 0)
		{
			ii = 1;
		}

		// Create the Row Matrix
		numat result(msize((ei - bi + 1) / ii, cols));
		for (int i = bi, k = 0; i <= ei; i += ii, k++)
		{
			for (int j = 0; j < cols; j++)
			{
				result(k, j) = (*this)(i, j);
			}
		}

		return result;
	}

	const numat numat::col(int idx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(idx >= 0 && idx < cols);

		// Get the selected column matrix
		numat result(msize(rows, 1));
		for (int i = 0; i < rows; i++)
		{
			result(i) = (*this)(i, idx);
		}

		return result;
	}

	const numat numat::col(int bi, int ei, int ii) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(bi >= 0 && bi < cols);
		assert(ei >= 0 && ei < cols);
		assert(bi < ei);
		if (ii < 0)
		{
			ii = 1;
		}

		// Create the Column Matrix
		numat result(msize(rows, (ei - bi + 1) / ii));
		for (int i = 0; i < rows; i++)
		{
			for (int j = bi, k = 0; j <= ei; j += ii, k++)
			{
				result(i, k) = (*this)(i, j);
			}
		}

		return result;
	}

	void numat::swapr(int idx0, int idx1) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the indices
		assert(idx0 >= 0 && idx0 < rows);
		assert(idx1 >= 0 && idx1 < rows);

		// Swap the matrix
		numat temp = row(idx0);
		for (int i = 0; i < cols; i++)
		{
			(*this)(idx0, i) = (*this)(idx1, i);
			(*this)(idx1, i) = temp(i);
		}
	}

	void numat::swapc(int idx0, int idx1) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the indices
		assert(idx0 >= 0 && idx0 < cols);
		assert(idx1 >= 0 && idx1 < cols);

		// Swap the matrix
		numat temp = col(idx0);
		for (int i = 0; i < rows; i++)
		{
			(*this)(i, idx0) = (*this)(i, idx1);
			(*this)(i, idx1) = temp(i);
		}
	}

	const numat numat::diag() const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the size
		assert(rows == cols);

		// Get a diagonal matrix
		numat result(msize(rows, cols), 0.0);
		for (int i = 0; i < rows; i++)
		{
			result(i, i) = (*this)(i, i);
		}

		return result;
	}

	numat numat::values(msize& size, double val)
	{
		// Create a values matrix
		return numat(size, val);
	}

	numat numat::zeros(msize& size)
	{
		// Create a zeros matrix
		return values(size, 0.0);
	}

	numat numat::ones(msize& size)
	{
		// Create an ones matrix
		return values(size, 1.0);
	}

	numat numat::eyes(int length)
	{
		// Create an eyes matrix
		numat result(msize(length, length), 0.0);
		for (int i = 0; i < length; i++)
		{
			result(i, i) = 1.0;
		}

		return result;
	}

	numat numat::permut(int length, int idx0, int idx1)
	{
		// Check the indices
		assert(idx0 >= 0 && idx0 < length);
		assert(idx1 >= 0 && idx1 < length);

		// Create a permutation matrix
		numat result = eyes(length);
		result(idx0, idx0) = 0.0;
		result(idx1, idx1) = 0.0;
		result(idx0, idx1) = 1.0;
		result(idx1, idx0) = 1.0;

		return result;
	}

	const numat numat::minorr(int idx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(idx >= 0 && idx < rows);

		// Check the size
		if (rows < 2)
		{
			return *this;
		}

		// Get a minor matrix
		numat result(msize(rows - 1, cols));
		for (int i = 0, k = 0; i < rows; i++)
		{
			// Check the index
			if (i == idx)
			{
				continue;
			}
			for (int j = 0; j < cols; j++)
			{
				result(k, j) = (*this)(i, j);
			}
			k++;
		}

		return result;
	}

	const numat numat::minorc(int idx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(idx >= 0 && idx < cols);

		// Check the size
		if (cols < 2)
		{
			return *this;
		}

		// Get a minor matrix
		numat result(msize(rows, cols - 1));
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0, k = 0; j < cols; j++)
			{
				// Check the index
				if (j == idx)
				{
					continue;
				}
				result(i, k) = (*this)(i, j);
				k++;
			}
		}

		return result;
	}

	const numat numat::minor(int ridx, int cidx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the indices
		assert(ridx >= 0 && ridx < rows);
		assert(cidx >= 0 && cidx < cols);

		// Get a minor matrix
		if (rows < 2 || cols < 2)
		{
			return *this;
		}
		else
		{
			return minorr(ridx).minorc(cidx);
		}
	}

	const double numat::expansion(const numat& mat) const
	{
		// Check the size
		double result = 0.0;
		if (mat.rows == 2 && mat.cols == 2)
		{
			// Calculate a determinant
			result = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
		}
		else
		{
			// Check the zero data
			int val = 0;
			int idx = 0;
			for (int i = 0; i < mat.rows; i++)
			{
				int zeros = 0;
				for (int j = 0; j < mat.cols; j++)
				{
					if (almost(mat(i, j), 0.0) == true)
					{
						zeros++;
					}
				}
				if (zeros > val)
				{
					val = zeros;
					idx = i;
				}
			}

			// Calculate a determinant
			for (int i = 0; i < mat.cols; i++)
			{
				if (almost(mat(idx, i), 0.0) == false)
				{
					result += std::pow(-1.0, idx + i) * mat(idx, i) * expansion(mat.minor(idx, i));
				}
			}
		}

		return result;
	}

	const numat numat::uniq() const
	{
		// Check the matrix
		assert(empty() == false);

		// Find the unique values in the matrix
		vector<double> vec;
		for (int i = 0; i < rows * cols; i++)
		{
			bool redundant = false;
			for (int j = 0; j < (int)vec.size(); j++)
			{
				if ((*this)(i) == vec[j])
				{
					redundant = true;
					break;
				}
			}
			if (redundant == false)
			{
				vec.push_back((*this)(i));
			}
		}

		// Get a unique matrix
		return to_mat(vec);
	}

	int numat::count(double val) const
	{
		// Check the matrix
		assert(empty() == false);

		// Count the value
		int result = 0;
		for (int i = 0; i < rows * cols; i++)
		{
			if ((*this)(i) == val)
			{
				result++;
			}
		}

		return result;
	}

	bool numat::almost(double val, double tar, int round)
	{
		// Compare the value
		double power = std::pow(10.0, round);
		long long _val = (long long)(val * power + 0.5);
		long long _tar = (long long)(tar * power + 0.5);
		if (_val == _tar)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	void numat::sortr(bool descending) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the sorting method
		if (descending == true)
		{
			// Sort the vector
			for (int i = 0; i < rows; i++)
			{
				int begin = i * cols;
				int end = (i + 1) * cols;
				std::sort(ptr() + begin, ptr() + end, greater<double>());
			}
		}
		else
		{
			// Sort the vector
			for (int i = 0; i < rows; i++)
			{
				int begin = i * cols;
				int end = (i + 1) * cols;
				std::sort(ptr() + begin, ptr() + end);
			}
		}
	}

	void numat::sortc(bool descending) const
	{
		// Check the matrix
		assert(empty() == false);

		// Sort the matrix
		numat temp = trans();
		temp.sortr(descending);
		for (int i = 0; i < temp.rows; i++)
		{
			for (int j = 0; j < temp.cols; j++)
			{
				(*this)(j, i) = temp(i, j);
			}
		}
	}

	void numat::sort(bool descending) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the sorting method
		if (descending == true)
		{
			// Sort the vector
			std::sort(ptr(), ptr() + length(), greater<double>());
		}
		else
		{
			// Sort the vector
			std::sort(ptr(), ptr() + length());
		}
	}

	numat numat::minr(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Find the minimum values
		numat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			double val = mat(i, 0);
			for (int j = 1; j < mat.cols; j++)
			{
				if (mat(i, j) < val)
				{
					val = mat(i, j);
				}
			}
			result(i) = val;
		}

		return result;
	}

	numat numat::minc(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Find the minimum values
		numat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			double val = mat(0, i);
			for (int j = 1; j < mat.rows; j++)
			{
				if (mat(j, i) < val)
				{
					val = mat(j, i);
				}
			}
			result(i) = val;
		}

		return result;
	}

	double numat::min(const numat& mat)
	{
		// Find a minimum value
		return minc(minr(mat))(0);
	}

	numat numat::maxr(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Find the maximum values
		numat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			double val = mat(i, 0);
			for (int j = 1; j < mat.cols; j++)
			{
				if (mat(i, j) > val)
				{
					val = mat(i, j);
				}
			}
			result(i) = val;
		}

		return result;
	}

	numat numat::maxc(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Find the maximum values
		numat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			double val = mat(0, i);
			for (int j = 1; j < mat.rows; j++)
			{
				if (mat(j, i) > val)
				{
					val = mat(j, i);
				}
			}
			result(i) = val;
		}

		return result;
	}

	double numat::max(const numat& mat)
	{
		// Find a maximum value
		return maxc(maxr(mat))(0);
	}

	numat numat::sumr(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the sum values
		numat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			double val = 0.0;
			for (int j = 0; j < mat.cols; j++)
			{
				val += mat(i, j);
			}
			result(i) = val;
		}

		return result;
	}

	numat numat::sumc(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the sum values
		numat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			double val = 0.0;
			for (int j = 0; j < mat.rows; j++)
			{
				val += mat(j, i);
			}
			result(i) = val;
		}

		return result;
	}

	double numat::sum(const numat& mat)
	{
		// Calculate a sum value
		return sumc(sumr(mat))(0);
	}

	numat numat::meanr(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the sum values
		numat temp = sumr(mat);

		// Calculate the mean values
		numat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			result(i) = temp(i) / mat.cols;
		}

		return result;
	}

	numat numat::meanc(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the sum values
		numat temp = sumc(mat);

		// Calculate the mean values
		numat result(msize(1, mat.cols));
		for (int i = 0; i < result.cols; i++)
		{
			result(i) = temp(i) / mat.rows;
		}

		return result;
	}

	double numat::mean(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a sum value
		double temp = sum(mat);

		// Calculate a mean value
		double result = 0.0;
		result = temp / mat.length();

		return result;
	}

	numat numat::varr(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the mean values
		numat temp = meanr(mat);

		// Calculate the variance values
		numat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			double val = 0.0;
			for (int j = 0; j < mat.cols; j++)
			{
				val += (mat(i, j) - temp(i)) * (mat(i, j) - temp(i));
			}
			result(i) = val / mat.cols;
		}

		return result;
	}

	numat numat::varc(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the mean values
		numat temp = meanc(mat);

		// Calculate the variance values
		numat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			double val = 0.0;
			for (int j = 0; j < mat.rows; j++)
			{
				val += (mat(j, i) - temp(i)) * (mat(j, i) - temp(i));
			}
			result(i) = val / mat.rows;
		}

		return result;
	}

	double numat::var(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a mean value
		double temp = mean(mat);

		// Calculate a Variance Value
		double result = 0.0;
		for (int i = 0; i < mat.rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				double val = mat(i, j) - temp;
				result += val * val;
			}
		}
		result /= mat.length();

		return result;
	}

	numat numat::stdr(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the variance values
		numat temp = varr(mat);

		// Calculate the standard deviation values
		numat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			result(i) = std::sqrt(temp(i));
		}

		return result;
	}

	numat numat::stdc(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate the variance values
		numat temp = varc(mat);

		// Calculate the standard deviation values
		numat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			result(i) = std::sqrt(temp(i));
		}

		return result;
	}

	double numat::std(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a variance value
		double temp = var(mat);

		// Calculate a standard deviation value
		double result = 0.0;
		result = std::sqrt(temp);

		return result;
	}

	const numat numat::cov() const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate the mean values
		numat mu = meanr(*this);

		// Calculate a covariance matrix
		numat result(msize(rows, 1), 0.0);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				result(i) += ((*this)(i, j) - mu(i)) * ((*this)(i, j) - mu(i));
			}
			result(i) /= cols;
		}

		return result;
	}

	void numat::log() const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate a logarithm matrix
		for (int i = 0; i < rows * cols; i++)
		{
			(*this)(i) = std::log((*this)(i));
		}
	}

	void numat::log(double base) const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate a logarithm matrix
		for (int i = 0; i < rows * cols; i++)
		{
			(*this)(i) = std::log((*this)(i)) / std::log(base);
		}
	}

	numat numat::log(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a logarithm matrix
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = std::log(mat(i));
		}

		return result;
	}

	numat numat::log(const numat& mat, double base)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a logarithm matrix
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = std::log(mat(i)) / std::log(base);
		}

		return result;
	}

	void numat::exp() const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate an exponential matrix
		for (int i = 0; i < rows * cols; i++)
		{
			(*this)(i) = std::exp((*this)(i));
		}
	}

	numat numat::exp(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate an exponential matrix
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = std::exp(mat(i));
		}

		return result;
	}

	void numat::pow(double p) const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate a power matrix
		for (int i = 0; i < rows * cols; i++)
		{
			(*this)(i) = std::pow((*this)(i), p);
		}
	}

	numat numat::pow(const numat& mat, double p)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a power matrix
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = std::pow(mat(i), p);
		}

		return result;
	}

	void numat::sqrt() const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate a square root matrix
		for (int i = 0; i < rows * cols; i++)
		{
			(*this)(i) = std::sqrt((*this)(i));
		}
	}

	numat numat::sqrt(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate a square root matrix
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = std::sqrt(mat(i));
		}

		return result;
	}

	void numat::abs() const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate an absolute matrix
		for (int i = 0; i < rows * cols; i++)
		{
			(*this)(i) = std::abs((*this)(i));
		}
	}

	numat numat::abs(const numat& mat)
	{
		// Check the matrix
		assert(mat.empty() == false);

		// Calculate an absolute matrix
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = std::abs(mat(i));
		}

		return result;
	}

	const numat numat::add(const numat& mat) const
	{
		// Check the matrix
		assert(empty() == false);
		assert(mat.empty() == false);

		// Check the size
		if (mat.length() == 1)
		{
			// Calculate an addition matrix
			return add(mat(0));
		}
		else if (length() == 1)
		{
			// Calculate an addition matrix
			return _ptr[0] + mat;
		}
		else if (rows == mat.rows && cols != mat.cols)
		{
			// Check the column size
			if (cols == 1)
			{
				// Calculate an addition matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i) + mat(i, j);
					}
				}

				return result;
			}
			else if (mat.cols == 1)
			{
				// Calculate an addition matrix
				numat result(msize(mat.rows, cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						result(i, j) = (*this)(i, j) + mat(i);
					}
				}

				return result;
			}
		}
		else if (rows != mat.rows && cols == mat.cols)
		{
			// Check the row size
			if (rows == 1)
			{
				// Calculate an addition matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(j) + mat(i, j);
					}
				}

				return result;
			}
			else if (mat.rows == 1)
			{
				// Calculate an addition matrix
				numat result(msize(rows, mat.cols));
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i, j) + mat(j);
					}
				}

				return result;
			}
		}

		// Calculate an addition matrix
		assert(rows == mat.rows && cols == mat.cols);
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = (*this)(i) + mat(i);
		}

		return result;
	}

	const numat numat::add(double val) const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate an addition matrix
		numat result(msize(rows, cols));
		for (int i = 0; i < rows * cols; i++)
		{
			result(i) = (*this)(i) + val;
		}

		return result;
	}

	const numat numat::sub(const numat& mat) const
	{
		// Check the matrix
		assert(empty() == false);
		assert(mat.empty() == false);

		// Check the size
		if (mat.length() == 1)
		{
			// Calculate a subtraction matrix
			return sub(mat(0));
		}
		else if (length() == 1)
		{
			// Calculate a subtraction matrix
			return _ptr[0] - mat;
		}
		else if (rows == mat.rows && cols != mat.cols)
		{
			// Check the column size
			if (cols == 1)
			{
				// Calculate a subtraction matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i) - mat(i, j);
					}
				}

				return result;
			}
			else if (mat.cols == 1)
			{
				// Calculate a subtraction matrix
				numat result(msize(mat.rows, cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						result(i, j) = (*this)(i, j) - mat(i);
					}
				}

				return result;
			}
		}
		else if (rows != mat.rows && cols == mat.cols)
		{
			// Check the row size
			if (rows == 1)
			{
				// Calculate a subtraction matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(j) - mat(i, j);
					}
				}

				return result;
			}
			else if (mat.rows == 1)
			{
				// Calculate a subtraction matrix
				numat result(msize(rows, mat.cols));
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i, j) - mat(j);
					}
				}

				return result;
			}
		}

		// Calculate a subtraction matrix
		assert(rows == mat.rows && cols == mat.cols);
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = (*this)(i) - mat(i);
		}

		return result;
	}

	const numat numat::sub(double val) const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate a subtraction matrix
		numat result(msize(rows, cols));
		for (int i = 0; i < rows * cols; i++)
		{
			result(i) = (*this)(i) - val;
		}

		return result;
	}

	const numat numat::mul(const numat& mat) const
	{
		// Check the matrix
		assert(empty() == false);
		assert(mat.empty() == false);

		// Check the size
		if (mat.length() == 1)
		{
			// Calculate a multiplication matrix
			return mul(mat(0));
		}
		else if (length() == 1)
		{
			// Calculate a multiplication matrix
			return _ptr[0] * mat;
		}
		else if (rows == mat.rows && cols != mat.cols)
		{
			// Check the column size
			if (cols == 1)
			{
				// Calculate a multiplication matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i) * mat(i, j);
					}
				}

				return result;
			}
			else if (mat.cols == 1)
			{
				// Calculate a multiplication matrix
				numat result(msize(mat.rows, cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						result(i, j) = (*this)(i, j) * mat(i);
					}
				}

				return result;
			}
		}
		else if (rows != mat.rows && cols == mat.cols)
		{
			// Check the row size
			if (rows == 1)
			{
				// Calculate a multiplication matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(j) * mat(i, j);
					}
				}

				return result;
			}
			else if (mat.rows == 1)
			{
				// Calculate a multiplication matrix
				numat result(msize(rows, mat.cols));
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i, j) * mat(j);
					}
				}

				return result;
			}
		}

		// Calculate a multiplication matrix
		assert(rows == mat.rows && cols == mat.cols);
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = (*this)(i) * mat(i);
		}

		return result;
	}

	const numat numat::mul(double val) const
	{
		// Check the matrix
		assert(empty() == false);

		// Calculate a multiplication matrix
		numat result(msize(rows, cols));
		for (int i = 0; i < rows * cols; i++)
		{
			result(i) = (*this)(i) * val;
		}

		return result;
	}

	const numat numat::div(const numat& mat) const
	{
		// Check the matrix
		assert(empty() == false);
		assert(mat.empty() == false);

		// Check the values
		for (int i = 0; i < mat.length(); i++)
		{
			assert(mat(i) != 0.0);
		}

		// Check the size
		if (mat.length() == 1)
		{
			// Calculate a division matrix
			return div(mat(0));
		}
		else if (length() == 1)
		{
			// Calculate a division matrix
			return _ptr[0] / mat;
		}
		else if (rows == mat.rows && cols != mat.cols)
		{
			// Check the column size
			if (cols == 1)
			{
				// Calculate a division matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i) / mat(i, j);
					}
				}

				return result;
			}
			else if (mat.cols == 1)
			{
				// Calculate a division matrix
				numat result(msize(mat.rows, cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						result(i, j) = (*this)(i, j) / mat(i);
					}
				}

				return result;
			}
		}
		else if (rows != mat.rows && cols == mat.cols)
		{
			// Check the row size
			if (rows == 1)
			{
				// Calculate a division matrix
				numat result(msize(mat.rows, mat.cols));
				for (int i = 0; i < mat.rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(j) / mat(i, j);
					}
				}

				return result;
			}
			else if (mat.rows == 1)
			{
				// Calculate a division matrix
				numat result(msize(rows, mat.cols));
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < mat.cols; j++)
					{
						result(i, j) = (*this)(i, j) / mat(j);
					}
				}

				return result;
			}
		}

		// Calculate a division matrix
		assert(rows == mat.rows && cols == mat.cols);
		numat result(msize(mat.rows, mat.cols));
		for (int i = 0; i < mat.rows * mat.cols; i++)
		{
			result(i) = (*this)(i) / mat(i);
		}

		return result;
	}

	const numat numat::div(double val) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the value
		assert(val != 0.0);

		// Calculate a division matrix
		numat result(msize(rows, cols));
		for (int i = 0; i < rows * cols; i++)
		{
			result(i) = (*this)(i) / val;
		}

		return result;
	}

	const numat numat::trans() const
	{
		// Check the matrix
		assert(empty() == false);

		// Get a transpose matrix
		numat result(msize(cols, rows));
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				result(j, i) = (*this)(i, j);
			}
		}

		return result;
	}

	const numat numat::T() const
	{
		// Get a transpose matrix
		return trans();
	}

	const numat numat::dot(const numat& mat) const
	{
		// Check the matrix
		assert(empty() == false);
		assert(mat.empty() == false);

		// Check the size
		assert(cols == mat.rows);

		// Get a dot product matrix
		numat result(msize(rows, mat.cols));
		for (int i = 0; i < result.rows; i++)
		{
			for (int j = 0; j < result.cols; j++)
			{
				double val = 0.0;
				for (int k = 0; k < cols; k++)
				{
					val += (*this)(i, k) * mat(k, j);
				}
				result(i, j) = val;
			}
		}

		return result;
	}

	const numat numat::cof() const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the size
		assert(rows == cols);

		// Get a cofactor expansion matrix
		numat result(msize(rows, cols));
		for (int i = 0; i < result.rows; i++)
		{
			for (int j = 0; j < result.cols; j++)
			{
				result(i, j) = std::pow(-1.0, i + j) * expansion(minor(i, j));
			}
		}

		return result;
	}

	const numat numat::inv() const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the size
		assert(rows == cols);

		// Get an inverse matrix
		numat result(msize(rows, cols));
		if (result.rows == 2 && result.cols == 2)
		{
			// Check a determinant
			double det = (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
			assert(almost(det, 0.0) == false);

			// Calculate the inverse Matrix
			result(0, 0) = (*this)(1, 1) / det;
			result(1, 1) = (*this)(0, 0) / det;
			result(1, 0) = -(*this)(1, 0) / det;
			result(0, 1) = -(*this)(0, 1) / det;
		}
		else
		{
			// Get a cafactor matrix
			numat c = cof();

			// Check a determinant
			double det = 0.0;
			for (int i = 0; i < cols; i++)
			{
				det += c(i) * (*this)(0, i);
			}
			assert(almost(det, 0.0) == false);

			// Calculate the inverse Matrix
			result = c.trans().mul(1 / det);
		}

		return result;
	}

	const numat numat::I() const
	{
		// Get an inverse matrix
		return inv();
	}

	const double numat::det() const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the size
		if (rows == 2 && cols == 2)
		{
			return expansion(*this);
		}
		else
		{
			// Get a cofactor matrix
			numat c(msize(1, cols));
			for (int i = 0; i < cols; i++)
			{
				c(i) = std::pow(-1.0, i) * expansion(minor(0, i));
			}

			// Calculate a determinant
			double result = 0.0;
			for (int i = 0; i < cols; i++)
			{
				result += c(i) * (*this)(0, i); 
			}

			return result;
		}
	}
}