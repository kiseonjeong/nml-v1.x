#include "stdafx.h"
#include "Spline.h"

namespace nml
{
	void spline::checkEmpty(const numat& X, const numat& Y) const
	{
		// Check the matrices are empty or not
		assert(X.empty() == false && Y.empty() == false);
	}

	void spline::checkLength(const numat& X, const numat& Y) const
	{
		// Check the length of the matrices
		assert(X.length() > 1 && Y.length() > 1);
	}

	void spline::compareLength(const numat& X, const numat& Y) const
	{
		// Compare the length between the matrices
		assert(X.length() == Y.length());
	}

	const numat spline::coeff() const
	{
		// Get the estimated coefficients
		return C;
	}

	const int spline::getInterval(const double x) const
	{
		// Check the interval matrix
		assert(I.empty() == false);

		// Check the left external region
		if (x < I(0, 0))
		{
			return 0;
		}
		// Check the right external region
		else if (x >= I(numit - 1, 1))
		{
			return numit - 1;
		}
		// Check the internal regions
		else
		{
			int idx = 0;
			for (int i = 0; i < numit; i++)
			{
				if (x >= I(i, 0) && x < I(i, 1))
				{
					idx = i;
					break;
				}
			}

			return idx;
		}
	}

	spline1d::spline1d()
	{
		// Set an object
		setObject();
	}

	spline1d::spline1d(const numat& X, const numat& Y)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y);
	}

	spline1d::spline1d(const spline1d& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	spline1d::~spline1d()
	{
		// Clear the object
		clearObject();
	}

	spline1d& spline1d::operator=(const spline1d& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void spline1d::setObject()
	{
		// Set the parameters
		setType(*this);
		numpt = 0;
		numit = 0;

		// Set the memories
		X.release();
		Y.release();
		C.release();
	}

	void spline1d::copyObject(const object& obj)
	{
		// Do down casting
		const spline1d* _obj = static_cast<const spline1d*>(&obj);

		// Copy the parameters
		numpt = _obj->numpt;
		numit = _obj->numit;

		// Copy the memories
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void spline1d::clearObject()
	{
		// Clear the memories
		X.release();
		Y.release();
		C.release();
	}

	void spline1d::estimate(const numat& X, const numat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the input vectors
		this->X = X.reshape(X.length());
		this->Y = Y.reshape(Y.length());

		// Set an interval matrix
		numpt = X.length();
		numit = numpt - 1;
		I = numat::zeros(msize(numit, 2));
		for (int i = 0; i < numit; i++)
		{
			I(i, 0) = X(i);
			I(i, 1) = X(i + 1);
		}

		// Calculate the coefficients
		calculate();
	}

	void spline1d::calculate()
	{
		// Calculate a divided difference value on each interval
		C = numat::zeros(msize(numit, 2));
		for (int i = 0; i < numit; i++)
		{
			C(i, 0) = Y(i);
			C(i, 1) = slope(X(i), Y(i), X(i + 1), Y(i + 1));
		}
	}

	const double spline1d::slope(const double x0, const double y0, const double x1, const double y1) const
	{
		// Calculate the slope
		return (y1 - y0) / (x1 - x0);
	}

	const double spline1d::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Get an interval index
		int i = getInterval(x);

		// Interpolate a Value on the input X
		double y = C(i, 0) + C(i, 1) * (x - X(i));

		return y;
	}

	spline2d::spline2d()
	{
		// Set an object
		setObject();
	}

	spline2d::spline2d(const numat& X, const numat& Y)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y);
	}

	spline2d::spline2d(const spline2d& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	spline2d::~spline2d()
	{
		// Clear the object
		clearObject();
	}

	spline2d& spline2d::operator=(const spline2d& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void spline2d::setObject()
	{
		// Set the parameters
		setType(*this);
	}

	void spline2d::estimate(const numat& X, const numat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the vectors
		this->X = X.reshape(X.length());
		this->Y = Y.reshape(Y.length());

		// Set an interval matrix
		numpt = X.length();
		numit = numpt - 1;
		I = numat::zeros(msize(numit, 2));
		for (int i = 0; i < numit; i++)
		{
			I(i, 0) = X(i);
			I(i, 1) = X(i + 1);
		}

		// Calculate the coefficients
		calculate();
	}

	void spline2d::calculate()
	{
		// Calculate an interval width
		numat W(msize(numit, 1));
		for (int i = 0; i < numit; i++)
		{
			W(i) = X(i + 1) - X(i);
		}

		// Calculate the slope
		numat S(msize(numit, 1));
		for (int i = 0; i < numit; i++)
		{
			S(i) = (Y(i + 1) - Y(i)) / W(i);
		}

		// Set the begin point conditions
		C = numat::zeros(msize(numit, 3));
		C(0, 2) = 0.0;
		C(0, 1) = S(0);

		// Calculate the 1st coefficients
		for (int i = 0; i < numit; i++)
		{
			C(i, 0) = Y(i);
		}

		// Calculate the 2nd and 3rd coefficients
		for (int i = 1; i < numit; i++)
		{
			C(i, 1) = C(i - 1, 1) + 2.0 * C(i - 1, 2) * W(i - 1);
			C(i, 2) = (S(i) - C(i, 1)) / W(i);
		}
	}

	const double spline2d::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Get an interval index
		int i = getInterval(x);

		// Interpolate a Value on the input X
		double y = C(i, 0) + C(i, 1) * (x - X(i)) + C(i, 2) * pow(x - X(i), 2);

		return y;
	}

	cubic_params::cubic_params() : type(_type)
	{
		// Set an object
		setObject();
	}

	cubic_params::~cubic_params()
	{
		// Clear the object
		clearObject();
	}

	void cubic_params::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = CUBIC_UNKNOWN;
	}

	void cubic_params::copyObject(const object& obj)
	{
		// Do down casting
		const cubic_params* _obj = static_cast<const cubic_params*>(&obj);

		// Copy the parameters
		_type = _obj->_type;
	}

	void cubic_params::clearObject()
	{

	}

	natural::natural()
	{
		// Set an object
		setObject();
	}

	natural::natural(const natural& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	natural::~natural()
	{
		// Clear the object
		clearObject();
	}

	natural& natural::operator=(const natural& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void natural::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = CUBIC_NATURAL;
	}

	clamped::clamped()
	{
		// Set an object
		setObject();
	}

	clamped::clamped(const clamped& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	clamped::~clamped()
	{
		// Clear the object
		clearObject();
	}

	clamped& clamped::operator=(const clamped& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void clamped::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = CUBIC_CLAMPED;
	}

	notaknot::notaknot() : alpha(_alpha), beta(_beta)
	{
		// Set an object
		setObject();
	}

	notaknot::notaknot(const double alpha, const double beta) : alpha(_alpha), beta(_beta)
	{
		// Set an object
		setObject();

		// Set the parameters
		set(alpha, beta);
	}

	notaknot::notaknot(const notaknot& obj) : alpha(_alpha), beta(_beta)
	{
		// Set an object
		setObject();

		// Clone Object
		*this = obj;
	}

	notaknot::~notaknot()
	{
		// Clear the object
		clearObject();
	}

	notaknot& notaknot::operator=(const notaknot& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void notaknot::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = CUBIC_NOTAKNOT;
		_alpha = 0.0;
		_beta = 0.0;
	}

	void notaknot::copyObject(const object& obj)
	{
		// Do down casting
		const notaknot* _obj = static_cast<const notaknot*>(&obj);

		// Copy the parameters
		_type = _obj->_type;
		_alpha = _obj->_alpha;
		_beta = _obj->_beta;
	}

	void notaknot::clearObject()
	{

	}

	void notaknot::set(const double alpha, const double beta)
	{
		// Set the parameters
		_alpha = alpha;
		_beta = beta;
	}

	spline3d::spline3d()
	{
		// Set an object
		setObject();
	}

	spline3d::spline3d(const cubic_params& params)
	{
		// Set an object
		setObject();

		// Set a cubic spline condition
		condition(params);
	}

	spline3d::spline3d(const numat& X, const numat& Y, const cubic_params& params)
	{
		// Set an object
		setObject();

		// Set a cubic spline condition
		condition(params);

		// Estimate the coefficients
		estimate(X, Y);
	}

	spline3d::spline3d(const spline3d& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	spline3d::~spline3d()
	{
		// Clear the object
		clearObject();
	}

	spline3d& spline3d::operator=(const spline3d& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void spline3d::setObject()
	{
		// Set the parameters
		setType(*this);
		numpt = 0;
		numit = 0;

		// Set the memories
		params = nullptr;
		X.release();
		Y.release();
		C.release();
	}

	void spline3d::copyObject(const object& obj)
	{
		// Do downcasting
		const spline3d* _obj = static_cast<const spline3d*>(&obj);

		// Copy the parameters
		numpt = _obj->numpt;
		numit = _obj->numit;

		// Copy the memories
		createParams(*_obj->params);
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void spline3d::clearObject()
	{
		// Clear the memories
		if (params != nullptr)
		{
			delete params;
		}
		X.release();
		Y.release();
		C.release();
	}

	void spline3d::createParams(const cubic_params& params)
	{
		// Check the cubic type and create the parameters
		switch (params.type)
		{
		case CUBIC_NATURAL: this->params = new natural((natural&)params); break;
		case CUBIC_CLAMPED: this->params = new clamped((clamped&)params); break;
		case CUBIC_NOTAKNOT: this->params = new notaknot((notaknot&)params); break;
		default: this->params = nullptr; break;
		}
	}

	void spline3d::condition(const cubic_params& params)
	{
		// Release the old parameters
		if (this->params != nullptr)
		{
			delete this->params;
		}

		// Create the parameters
		createParams(params);
	}

	void spline3d::estimate(const numat& X, const numat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the input vectors
		this->X = X.reshape(X.length());
		this->Y = Y.reshape(Y.length());

		// Set an interval matrix
		numpt = X.length();
		numit = numpt - 1;
		I = numat::zeros(msize(numit, 2));
		for (int i = 0; i < numit; i++)
		{
			I(i, 0) = X(i);
			I(i, 1) = X(i + 1);
		}

		// Calculate the coefficients
		calculate();
	}

	void spline3d::calculate()
	{
		// Check the cubic type
		assert(params->type != CUBIC_UNKNOWN);

		// Calculate an interval width
		numat W(msize(numit, 1));
		for (int i = 0; i < numit; i++)
		{
			W(i) = X(i + 1) - X(i);
		}

		// Calculate the slope
		numat S(msize(numit, 1));
		for (int i = 0; i < numit; i++)
		{
			S(i) = (Y(i + 1) - Y(i)) / W(i);
		}

		// Get a condition matrix
		numat cmat = condmat(W);

		// Get a result matrix
		numat rmat = resmat(S);

		// Solve the condition matrix
		numat smat = cmat.I() * rmat;

		// Calculate the 1st coefficients
		C = numat::zeros(msize(numit, 4));
		for (int i = 0; i < numit; i++)
		{
			C(i, 0) = Y(i);
		}

		// Calculate the 2nd coefficients
		for (int i = 0; i < numit; i++)
		{
			C(i, 1) = S(i) - W(i) / 3.0 * (2.0 * smat(i) + smat(i + 1));
		}

		// Calculate the 3rd coefficients
		for (int i = 0; i < numit; i++)
		{
			C(i, 2) = smat(i);
		}

		// Calculate the 4th coefficients
		for (int i = 0; i < numit; i++)
		{
			C(i, 3) = (smat(i + 1) - smat(i)) / (3.0 * W(i));
		}
	}

	const numat spline3d::condmat(const numat& W) const
	{
		// Create a condition matrix
		numat cond(msize(numpt, numpt), 0.0);
		for (int i = 1; i < cond.rows - 1; i++)
		{
			cond(i, i - 1) = W(i - 1);
			cond(i, i) = 2.0 * (W(i - 1) + W(i));
			cond(i, i + 1) = W(i);
		}

		// Apply the end point conditions
		switch (params->type)
		{
		case CUBIC_NATURAL:
			cond(0, 0) = 1.0;
			cond(0, 1) = 0.0;
			cond(0, 2) = 0.0;
			cond(numpt - 1, numpt - 3) = 0.0;
			cond(numpt - 1, numpt - 2) = 0.0;
			cond(numpt - 1, numpt - 1) = 1.0;
			break;
		case CUBIC_CLAMPED:
			cond(0, 0) = 2.0 * W(0);
			cond(0, 1) = W(0);
			cond(0, 2) = 0.0;
			cond(numpt - 1, numpt - 3) = 0.0;
			cond(numpt - 1, numpt - 2) = W(numit - 1);
			cond(numpt - 1, numpt - 1) = 2.0 * W(numit - 1);
			break;
		case CUBIC_NOTAKNOT:
			cond(0, 0) = -W(1);
			cond(0, 1) = W(0) + W(1);
			cond(0, 2) = -W(0);
			cond(numpt - 1, numpt - 3) = -W(numit - 1);
			cond(numpt - 1, numpt - 2) = W(numit - 1) + W(numit - 2);
			cond(numpt - 1, numpt - 1) = -W(numit - 2);
			break;
		default:
			cond(0, 0) = 1.0;
			cond(0, 1) = 0.0;
			cond(0, 2) = 0.0;
			cond(numpt - 1, numpt - 3) = 0.0;
			cond(numpt - 1, numpt - 2) = 0.0;
			cond(numpt - 1, numpt - 1) = 1.0;
			break;
		}

		return cond;
	}

	const numat spline3d::resmat(const numat& S) const
	{
		// Create a result matrix
		numat res(msize(numpt, 1), 0.0);
		for (int i = 1; i < res.length() - 1; i++)
		{
			res(i) = 3.0 * (S(i) - S(i - 1));
		}

		// Apply the end point conditions
		switch (params->type)
		{
		case CUBIC_NATURAL:
		case CUBIC_CLAMPED:
			res(0) = 0.0;
			res(numpt - 1) = 0.0;
			break;
		case CUBIC_NOTAKNOT:
			res(0) = S(0) - ((notaknot*)params)->alpha;
			res(numpt - 1) = ((notaknot*)params)->beta - S(numit - 1);
			break;
		default:
			res(0) = 0.0;
			res(numpt - 1) = 0.0;
			break;
		}

		return res;
	}

	const double spline3d::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Get an interval index
		int i = getInterval(x);

		// Interpolate a value on the input X
		double y = C(i, 0) + C(i, 1) * (x - X(i)) + C(i, 2) * pow(x - X(i), 2) + C(i, 3) * pow(x - X(i), 3);

		return y;
	}

	hermite_params::hermite_params() : type(_type)
	{
		// Set an object
		setObject();
	}

	hermite_params::~hermite_params()
	{
		// Clear the object
		clearObject();
	}

	void hermite_params::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = HERMITE_UNKNOWN;
	}

	void hermite_params::copyObject(const object& obj)
	{
		// Do down casting
		const hermite_params* _obj = static_cast<const hermite_params*>(&obj);

		// Copy the parameters
		_type = _obj->_type;
	}

	void hermite_params::clearObject()
	{

	}

	finitediff::finitediff()
	{
		// Set an object
		setObject();
	}

	finitediff::finitediff(const finitediff& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	finitediff::~finitediff()
	{
		// Clear the object
		clearObject();
	}

	finitediff& finitediff::operator=(const finitediff& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void finitediff::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = HERMITE_FINITE_DIFFERENCE;
	}

	const numat finitediff::calculate(const numat& X, const numat& Y, const int numpt) const
	{
		// Calculate the begin point condition
		numat M(msize(numpt, 1), 0.0);
		M(0) = slope(X(0), Y(0), X(1), Y(1));

		// Calculate the end point condition
		M(numpt - 1) = slope(X(numpt - 2), Y(numpt - 2), X(numpt - 1), Y(numpt - 1));

		// Calculate the others
		for (int i = 1; i < numpt - 1; i++)
		{
			M(i) = (slope(X(i), Y(i), X(i + 1), Y(i + 1)) + slope(X(i), Y(i), X(i - 1), Y(i - 1))) / 2.0;
		}

		return M;
	}

	const double finitediff::slope(const double x0, const double y0, const double x1, const double y1) const
	{
		// Calculate the slope
		return (y1 - y0) / (x1 - x0);
	}

	cardinal::cardinal() : c(_c)
	{
		// Set an object
		setObject();
	}

	cardinal::cardinal(double c) : c(_c)
	{
		// Set an object
		setObject();

		// Set the parameters
		set(c);
	}

	cardinal::cardinal(const cardinal& obj) : c(_c)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	cardinal::~cardinal()
	{
		// Clear the object
		clearObject();
	}

	cardinal& cardinal::operator=(const cardinal& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void cardinal::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = HERMITE_CARDINAL;
		_c = 0.0;
	}

	void cardinal::copyObject(const object& obj)
	{
		// Do down casting
		const cardinal* _obj = static_cast<const cardinal*>(&obj);

		// Copy the parameters
		_type = _obj->_type;
		_c = _obj->_c;
	}

	void cardinal::set(double c)
	{
		// Check the tension parameter
		assert(c >= 0.0 && c <= 1.0);

		// Set the parameter
		_c = c;
	}

	const numat cardinal::calculate(const numat& X, const numat& Y, const int numpt) const
	{
		// Calculate the begin point condition
		numat M(msize(numpt, 1), 0.0);
		M(0) = (1 - _c) * slope(X(0), Y(0), X(1), Y(1));

		// Calculate the end point condition
		M(numpt - 1) = (1 - _c) * slope(X(numpt - 2), Y(numpt - 2), X(numpt - 1), Y(numpt - 1));

		// Calculate the others
		for (int i = 1; i < numpt - 1; i++)
		{
			M(i) = (1 - _c) * slope(X(i - 1), Y(i - 1), X(i + 1), Y(i + 1));
		}

		return M;
	}

	catmullrom::catmullrom()
	{
		// Set an object
		setObject();
	}

	catmullrom::catmullrom(const catmullrom& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	catmullrom::~catmullrom()
	{
		// Clear the object
		clearObject();
	}

	catmullrom& catmullrom::operator=(const catmullrom& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void catmullrom::setObject()
	{
		// Set the parameters
		setType(*this);
		_type = HERMITE_CATMULL_ROM;
	}

	const numat catmullrom::calculate(const numat& X, const numat& Y, const int numpt) const
	{
		// Calculate the begin point condition
		numat M(msize(numpt, 1), 0.0);
		M(0) = (1 - _c) * slope(X(0), Y(0), X(1), Y(1));

		// Calculate the end point condition
		M(numpt - 1) = (1 - _c) * slope(X(numpt - 2), Y(numpt - 2), X(numpt - 1), Y(numpt - 1));

		// Calculate the others
		for (int i = 1; i < numpt - 1; i++)
		{
			M(i) = (1 - _c) * slope(X(i - 1), Y(i - 1), X(i + 1), Y(i + 1));
		}

		return M;
	}

	hermite::hermite()
	{
		// Set an object
		setObject();
	}

	hermite::hermite(const hermite_params& params)
	{
		// Set an object
		setObject();

		// Set a hermite spline condition
		condition(params);
	}

	hermite::hermite(const numat& X, const numat& Y, const hermite_params& params)
	{
		// Set an object
		setObject();

		// Set a hermite spline condition
		condition(params);

		// Estimate the coefficients
		estimate(X, Y);
	}

	hermite::hermite(const hermite& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	hermite::~hermite()
	{
		// Clear the object
		clearObject();
	}

	hermite& hermite::operator=(const hermite& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void hermite::setObject()
	{
		// Set the parameters
		setType(*this);
		numpt = 0;
		numit = 0;

		// Set the memories
		params = nullptr;
		X.release();
		Y.release();
		C.release();
	}

	void hermite::copyObject(const object& obj)
	{
		// Do down casting
		const hermite* _obj = static_cast<const hermite*>(&obj);

		// Copy the parameters
		numpt = _obj->numpt;
		numit = _obj->numit;

		// Copy the memories
		createParams(*_obj->params);
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void hermite::clearObject()
	{
		// Clear the memories
		if (params != nullptr)
		{
			delete params;
		}
		X.release();
		Y.release();
		C.release();
	}

	void hermite::createParams(const hermite_params& params)
	{
		// Check the hermite type and create the parameters
		switch (params.type)
		{
		case HERMITE_FINITE_DIFFERENCE: this->params = new finitediff((finitediff&)params); break;
		case HERMITE_CARDINAL: this->params = new cardinal((cardinal&)params); break;
		case HERMITE_CATMULL_ROM: this->params = new catmullrom((catmullrom&)params); break;
		default: this->params = nullptr; break;
		}
	}

	void hermite::condition(const hermite_params& params)
	{
		// Release the old parameters
		if (this->params != nullptr)
		{
			delete this->params;
		}

		// Create the parameters
		createParams(params);
	}

	void hermite::estimate(const numat& X, const numat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the iuput vectors
		this->X = X.reshape(X.length());
		this->Y = Y.reshape(Y.length());

		// Set an interval matrix
		numpt = X.length();
		numit = numpt - 1;
		I = numat::zeros(msize(numit, 2));
		for (int i = 0; i < numit; i++)
		{
			I(i, 0) = X(i);
			I(i, 1) = X(i + 1);
		}

		// Check the hermite spline condition
		assert(params != nullptr);

		// Calculate the coefficients
		C = params->calculate(this->X, this->Y, numpt);
	}

	const double hermite::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Get an interval index
		int i = getInterval(x);

		// Calculate a hermite basis values
		double h00 = basis(x, i, 0);
		double h10 = basis(x, i, 1);
		double h01 = basis(x, i, 2);
		double h11 = basis(x, i, 3);

		// Interpolate a value on the input X
		double y = h00 * Y(i) + h10 * C(i) + h01 * Y(i + 1) + h11 * C(i + 1);

		return y;
	}

	const double hermite::basis(const double x, const int idx, const int mode) const
	{
		// Calculate a tangent Value
		double t = tangent(x, idx);

		// Calculate a basis value
		double h = 0.0;
		switch (mode)
		{
		case 0: h = 2.0 * pow(t, 3.0) - 3.0 * pow(t, 2.0) + 1; break;
		case 1: h = pow(t, 3.0) - 2.0 * pow(t, 2.0) + t; break;
		case 2: h = -2.0 * pow(t, 3.0) + 3.0 * pow(t, 2.0); break;
		case 3: h = pow(t, 3.0) - pow(t, 2.0); break;
		}

		return h;
	}

	const double hermite::tangent(const double x, const int idx) const
	{
		// Calculate a tangent for an arbitary interval
		return (x - X(idx)) / (X(idx + 1) - X(idx));
	}

	akima::akima()
	{
		// Set an object
		setObject();
	}

	akima::akima(const numat& X, const numat& Y)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y);
	}

	akima::akima(const akima& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	akima::~akima()
	{
		// Clear the object
		clearObject();
	}

	akima& akima::operator=(const akima& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void akima::setObject()
	{
		// Set the parameters
		setType(*this);
		numpt = 0;
		numit = 0;

		// Clear the object
		clearObject();
	}

	void akima::copyObject(const object& obj)
	{
		// Do down casting
		const akima* _obj = static_cast<const akima*>(&obj);

		// Copy the parameters
		numpt = _obj->numpt;
		numit = _obj->numit;

		// Copy the memories
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void akima::clearObject()
	{
		// Clear the memories
		X.release();
		Y.release();
		C.release();
	}

	void akima::estimate(const numat& X, const numat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the input vectors
		this->X = X.reshape(X.length());
		this->Y = Y.reshape(Y.length());

		// Set an interval matrix
		numpt = X.length();
		numit = numpt - 1;
		I = numat::zeros(msize(numit, 2));
		for (int i = 0; i < numit; i++)
		{
			I(i, 0) = X(i);
			I(i, 1) = X(i + 1);
		}

		// Calculate the coefficients
		calculate();
	}

	void akima::calculate()
	{
		// Calculate the slopes on each interval
		numat dX(msize(numit, 1));
		numat dY(msize(numit, 1));
		numat S(msize(numit, 1));
		for (int i = 0; i < numit; i++)
		{
			dX(i) = X(i + 1) - X(i);
			dY(i) = Y(i + 1) - Y(i);
			S(i) = dY(i) / dX(i);
		}

		// Augment the slopes
		numat A = augment(S);

		// Calculate the coefficients
		C = numat::zeros(msize(numit, 4));
		for (int i = 0; i < numit; i++)
		{
			double dL = derivative(A, i);
			double dR = derivative(A, i + 1);
			C(i, 0) = Y(i);
			C(i, 1) = dL;
			C(i, 2) = (3.0 * A(i + 2) - dR - 2.0 * dL) / dX(i);
			C(i, 3) = -(2.0 * A(i + 2) - dR - dL) / pow(dX(i), 2.0);
		}
	}

	const numat akima::augment(const numat& S) const
	{
		// Calculate the begin point condition
		numat B(msize(2, 1));
		B(0) = 2.0 * (2.0 * S(0) - S(1)) - S(0);
		B(1) = 2.0 * S(0) - S(1);

		// Calculate the end point condition
		numat E(msize(2, 1));
		E(0) = 2.0 * S(numit - 1) - S(numit - 2);
		E(1) = 2.0 * (2.0 * S(numit - 1) - S(numit - 2)) - S(numit - 1);

		// Augment the vectors
		return B.vappend(S).vappend(E);
	}

	const double akima::derivative(const numat& S, int idx) const
	{
		// Calculate a derivative value
		double dL = abs(S(idx + 3) - S(idx + 2));
		double dR = abs(S(idx + 1) - S(idx));
		double dS = dL + dR;
		double wL = dL / dS;
		double wR = dR / dS;
		if (dL == 0.0 && dR == 0.0)
		{
			return (S(idx + 1) + S(idx + 2)) / 2.0;
		}
		else
		{
			return wL * S(idx + 1) + wR * S(idx + 2);
		}
	}

	const double akima::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Get an interval index
		int i = getInterval(x);

		// Interpolate a value on the input X
		double y = C(i, 0) + C(i, 1) * (x - X(i)) + C(i, 2) * pow(x - X(i), 2) + C(i, 3) * pow(x - X(i), 3);

		return y;
	}
}