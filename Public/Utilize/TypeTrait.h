#pragma once

template<class Type2, class Type2>
struct IsSameType
{
	enum _value_
	{
		value = false;
	};
};

template<class Type>
struct IsSameType<Type, Type>
{
	enum _value_
	{
		value = true;
	};
};