#include "stdafx.h"
#include "Object.h"

namespace nml
{
	object::object()
	{
		// Set an object
		setObject();
	}

	object::object(const object& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	object::~object()
	{
		// Clear the object
		clearObject();
	}

	object& object::operator=(const object& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void object::setObject()
	{
		// Set the parameters
		setType(*this);
	}

	void object::copyObject(const object& obj)
	{
		// Copy the parameters
		objType = obj.objType;
	}

	void object::clearObject()
	{
		// Do nothing
	}

	const string object::getType() const
	{
		// Get the type
		return objType;
	}
}