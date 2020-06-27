#include "stdafx.h"
#include "Regression.h"

namespace nml
{
	void regression::checkEmpty(const numat& X, const numat& Y) const
	{
		// Check the matrices are empty or not
		assert(X.empty() == false && Y.empty() == false);
	}

	void regression::checkLength(const numat& X, const numat& Y) const
	{
		// Check the length of the matrices
		assert(X.length() > 1 && Y.length() > 1);
	}

	void regression::compareLength(const numat& X, const numat& Y) const
	{
		// Compare the length between the matrices
		assert(X.length() == Y.length());
	}

	const numat regression::coeff() const
	{
		// Get the estimated coefficients
		return C;
	}

	curvefit::curvefit()
	{
		// Set an object
		setObject();
	}

	curvefit::curvefit(const numat& X, const numat& Y, int deg)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y, deg);
	}

	curvefit::curvefit(const curvefit& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	curvefit::~curvefit()
	{
		// Clear the object
		clearObject();
	}

	curvefit& curvefit::operator=(const curvefit& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void curvefit::setObject()
	{
		// Set the parameters
		setType(*this);
		epsilon = 0.01;
		meanX = 0.0;
		stdX = 1.0;
		meanY = 0.0;
		stdY = 1.0;

		// Set the memories
		C.release();
	}

	void curvefit::copyObject(const object& obj)
	{
		// Do down casting
		curvefit* _obj = (curvefit*)&obj;

		// Copy the parameters
		deg = _obj->deg;
		epsilon = _obj->epsilon;
		meanX = _obj->meanX;
		stdX = _obj->stdX;
		meanY = _obj->meanY;
		stdY = _obj->stdY;

		// Copy the memories
		C = _obj->C;
	}

	void curvefit::clearObject()
	{
		// Clear the memories
		C.release();
	}

	void curvefit::estimate(const numat& X, const numat& Y, const int deg)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		C.release();

		// Backup the input vectors
		this->X = X.reshape(X.length());
		this->Y = Y.reshape(Y.length());

		// Set the degree
		this->deg = deg;

		// Do whitening on the vectors
		whitening();

		// Calculate the coefficients
		calculate();
	}

	void curvefit::whitening()
	{
		// Do whitening the vector X
		meanX = numat::mean(X);
		stdX = numat::std(X);
		X = (X - meanX) / (stdX + epsilon);

		// Do whitening the vector Y
		meanY = numat::mean(Y);
		stdY = numat::std(Y);
		Y = (Y - meanY) / (stdY + epsilon);
	}

	void curvefit::calculate()
	{
		// Create a normal equation on the vector X
		numat Xp(msize(deg + 1, deg + 1), 0.0);
		for (int i = 0; i < Xp.rows; i++)
		{
			for (int j = 0; j < Xp.cols; j++)
			{
				for (int k = 0; k < X.length(); k++)
				{
					Xp(i, j) += pow(X(k), i + j);
				}
			}
		}

		// Create a normal equation on the vector Y
		numat Yp(msize(deg + 1, 1), 0.0);
		for (int i = 0; i < Yp.rows; i++)
		{
			for (int j = 0; j < Y.length(); j++)
			{
				Yp(i) += pow(X(j), i) * Y(j);
			}
		}

		// Calculate the coefficients
		C = SVD::pinv(Xp) * Yp;
	}

	const double curvefit::fit(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Fit a value on the input X
		double xp = (x - meanX) / (stdX + epsilon);
		double y = 0.0;
		for (int i = 0; i < C.length(); i++)
		{
			y += C(i) * pow(xp, i);
		}

		return y * (stdY + epsilon) + meanY;
	}
}