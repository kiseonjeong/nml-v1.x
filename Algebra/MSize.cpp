#include "stdafx.h"
#include "MSize.h"

namespace nml
{
	msize::msize()
	{
		d0 = -1;
		d1 = -1;
		d2 = -1;
	}

	msize::msize(const int length)
	{
		// Set a 1D size information
		set(length);
	}

	msize::msize(const int rows, const int cols)
	{
		// Set a 2D size information
		set(rows, cols);
	}

	msize::msize(const int rows, const int cols, const int channels)
	{
		// Set a 3D size information
		set(rows, cols, channels);
	}

	msize::~msize()
	{

	}
	
	msize& msize::operator=(const msize& obj)
	{
		// Copy from the input object
		d0 = obj.d0;
		d1 = obj.d1;
		d2 = obj.d2;

		return *this;
	}

	void msize::set(const int length)
	{
		// Set the dimension parameters
		d0 = length;
		d1 = 1;
		d2 = 1;
	}

	void msize::set(const int rows, const int cols)
	{
		// Set the dimension parameters
		d0 = rows;
		d1 = cols;
		d2 = 1;
	}

	void msize::set(const int rows, const int cols, const int channels)
	{
		// Set the dimension parameters
		d0 = rows;
		d1 = cols;
		d2 = channels;
	}
}