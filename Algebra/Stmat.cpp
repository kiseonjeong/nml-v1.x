#include "stdafx.h"
#include "Stmat.h"

namespace nml
{
	stmat::stmat()
	{
		// Set an object
		setObject();
	}

	stmat::stmat(const msize& size)
	{
		// Set an object
		setObject();

		// Create a matrix
		create(size);
	}

	stmat::stmat(const msize& size, string val)
	{
		// Set an object
		setObject();

		// Create a matrix
		create(size, val);
	}

	stmat::stmat(const numem<string>& data)
	{
		// Set an object
		setObject();

		// Create a matrix
		copy(data);
	}

	stmat::stmat(const stmat& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	stmat::~stmat()
	{
		// Clear the object
		clearObject();
	}

	stmat& stmat::operator=(const stmat& mat)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(mat);

		return *this;
	}

	void stmat::setObject()
	{
		// Set the parameters
		setType(*this);
	}

	void stmat::cout() const
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

	const numem<string> stmat::to_data() const
	{
		// Check the matrix
		if (empty() == true)
		{
			return numem<string>();
		}

		// Convert the matrix to the 'numem' data
		return numem<string>(*this);
	}

	stmat stmat::to_mat(const numem<string>& data)
	{
		// Check the data
		if (data.empty() == true)
		{
			return stmat();
		}

		// Convert the data to the matrix
		return stmat(data);
	}

	stmat stmat::to_mat(const vector<string>& vec)
	{
		// Check the vector
		if (vec.empty() == true)
		{
			return stmat();
		}

		// Create a matrix
		stmat result(msize((int)vec.size(), 1));

		// Initialize the matrix
		for (int i = 0; i < result.length(); i++)
		{
			result(i) = vec[i];
		}

		return result;
	}

	stmat stmat::to_mat(const vector<string>& vec, int nrows)
	{
		// Check the vector
		if (vec.empty() == true)
		{
			return stmat();
		}

		// Convert the vector to the matrix
		stmat temp = to_mat(vec);

		// Reshape the matrix
		return temp.reshape(nrows);
	}

	const stmat stmat::clone() const
	{
		// Check the matrix
		if (empty() == true)
		{
			return stmat();
		}

		// Clone the matrix
		return stmat(*this);
	}

	const stmat stmat::reshape(int nrows) const
	{
		// Check the matrix
		if (empty() == true)
		{
			return stmat();
		}

		// Check the size
		assert(nrows > 0 && nrows <= length());

		// Check the dimension
		assert(length() % nrows == 0);

		// Check the current row information
		if (rows == nrows)
		{
			return stmat(*this);
		}

		// Set the rows and the columns
		int rows = nrows;
		int cols = (int)ceil((double)length() / nrows);

		// Create a reshaped matrix
		stmat result(msize(rows, cols), "");
		for (int i = 0; i < length(); i++)
		{
			result(i) = (*this)(i);
		}

		return result;
	}

	stmat stmat::vappend(const vector<stmat>& arr)
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
		stmat result(msize(rows, cols));
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

	const stmat stmat::vappend(const stmat& mat) const
	{
		// Create an array
		vector<stmat> arr;
		arr.push_back(*this);
		arr.push_back(mat);

		// Append the matrices
		return vappend(arr);
	}

	stmat stmat::happend(const vector<stmat>& arr)
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
		stmat result(msize(rows, cols));
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

	const stmat stmat::happend(const stmat& mat) const
	{
		// Create an array
		vector<stmat> arr;
		arr.push_back(*this);
		arr.push_back(mat);

		// Append the matrices
		return happend(arr);
	}

	const stmat stmat::row(int idx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(idx >= 0 && idx < rows);

		// Get the selected row matrix
		stmat result(msize(1, cols));
		for (int i = 0; i < cols; i++)
		{
			result(i) = (*this)(idx, i);
		}

		return result;
	}

	const stmat stmat::col(int idx) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the index
		assert(idx >= 0 && idx < cols);

		// Get the selected column matrix
		stmat result(msize(rows, 1));
		for (int i = 0; i < rows; i++)
		{
			result(i) = (*this)(i, idx);
		}

		return result;
	}

	void stmat::swapr(int idx0, int idx1) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the indices
		assert(idx0 >= 0 && idx0 < rows);
		assert(idx1 >= 0 && idx1 < rows);

		// Swap the matrix
		stmat temp(*this);
		for (int i = 0; i < cols; i++)
		{
			(*this)(idx0, i) = temp(idx1, i);
			(*this)(idx1, i) = temp(idx0, i);
		}
	}

	void stmat::swapc(int idx0, int idx1) const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the indices
		assert(idx0 >= 0 && idx0 < cols);
		assert(idx1 >= 0 && idx1 < cols);

		// Swap the matrix
		stmat temp(*this);
		for (int i = 0; i < rows; i++)
		{
			(*this)(i, idx0) = temp(i, idx1);
			(*this)(i, idx1) = temp(i, idx0);
		}
	}

	const stmat stmat::diag() const
	{
		// Check the matrix
		assert(empty() == false);

		// Check the size
		assert(rows == cols);

		// Get a diagonal matrix
		stmat result(msize(rows, cols), "");
		for (int i = 0; i < rows; i++)
		{
			result(i, i) = (*this)(i, i);
		}

		return result;
	}

	const stmat stmat::minorr(int idx) const
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
		stmat result(msize(rows - 1, cols));
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

	const stmat stmat::minorc(int idx) const
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
		stmat result(msize(rows, cols - 1));
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

	const stmat stmat::minor(int ridx, int cidx) const
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

	const stmat stmat::uniq() const
	{
		// Check the matrix
		assert(empty() == false);

		// Find the unique values in the matrix
		vector<string> vec;
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

	int stmat::count(string val) const
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
}