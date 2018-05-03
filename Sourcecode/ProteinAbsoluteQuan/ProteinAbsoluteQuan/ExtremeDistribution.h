/*
#  Copyright(C) 2015-2018 all rights reserved

#  This program is a free software; you can redistribute it and / or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/Licenses/
*/

#pragma once
#include<random>

#pragma pack(push,_CRT_PACKING)
#pragma warning(push,3)
#pragma push_macro("new")
#undef new

#pragma warning(disable: 4127 4244 4296 4389 4521 4723)

#pragma warning(disable: 6294 6326)

#ifndef _BITS_BYTE
#define _BITS_BYTE	8
#endif /* _BITS_BYTE */

#if defined(_DEBUG) || defined(_RNG_CHECK)
#define _RNG_ASSERT(ex, msg)	\
	((ex) ? (void)0 : _Rng_abort(__FILE__ "(" _STRINGIZE(__LINE__) "): " msg))

#else /* defined(_DEBUG) || defined(_RNG_CHECK) */
#define _RNG_ASSERT(ex, msg) ((void)0)
#endif /* defined(_DEBUG) || defined(_RNG_CHECK) */

_STD_BEGIN



template<class _Ty = double>
class extreme_value_distribution_MinCase // Change By GZHQ 20170331
{	// template class for extreme value distribution
public:
	//static_assert(_Is_RealType<_Ty>::value,
	//	"invalid template argument for extreme_value_distribution");

	typedef extreme_value_distribution_MinCase<_Ty> _Myt;
	typedef _Ty result_type;

	struct param_type
	{	// parameter package
		typedef _Myt distribution_type;

		param_type(_Ty _A0 = _Ty(0),
			_Ty _B0 = _Ty(1))
		{	// construct from parameters
			_Init(_A0, _B0);
		}

		bool operator==(const param_type& _Right) const
		{	// test for equality
			return (_Ax == _Right._Ax && _Bx == _Right._Bx);
		}

		bool operator!=(const param_type& _Right) const
		{	// test for inequality
			return (!(*this == _Right));
		}

		_Ty a() const
		{	// return a value
			return (_Ax);
		}

		_Ty b() const
		{	// return b value
			return (_Bx);
		}

		void _Init(_Ty _A0, _Ty _B0)
		{	// initialize
			_RNG_ASSERT(0.0 < _B0,
				"invalid b argument for extreme_value_distribution");
			_Ax = _A0;
			_Bx = _B0;
		}

		_Ty _Ax;
		_Ty _Bx;
	};

	explicit extreme_value_distribution_MinCase(_Ty _A0 = _Ty(0),
		_Ty _B0 = _Ty(1))
		: _Par(_A0, _B0)
	{	// construct
	}

	explicit extreme_value_distribution_MinCase(param_type _Par0)
		: _Par(_Par0)
	{	// construct from parameter package
	}

	_Ty a() const
	{	// return a value
		return (_Par.a());
	}

	_Ty b() const
	{	// return b value
		return (_Par.b());
	}

	param_type param() const
	{	// return parameter package
		return (_Par);
	}

	void param(const param_type& _Par0)
	{	// set parameter package
		_Par = _Par0;
	}

	result_type(min)() const
	{	// get smallest possible result
		return (-(numeric_limits<result_type>::max)());
	}

	result_type(max)() const
	{	// get largest possible result
		return ((numeric_limits<result_type>::max)());
	}

	void reset()
	{	// clear internal state
	}

	template<class _Engine>
	result_type operator()(_Engine& _Eng) const
	{	// return next value
		return (_Eval(_Eng, _Par));
	}

	template<class _Engine>
	result_type operator()(_Engine& _Eng, const param_type& _Par0) const
	{	// return next value, given parameter package
		return (_Eval(_Eng, _Par0));
	}

	template<class _Elem,
	class _Traits>
		basic_istream<_Elem, _Traits>& _Read(
		basic_istream<_Elem, _Traits>& _Istr)
	{	// read state from _Istr
		_Ty _A0;
		_Ty _B0;
		_In(_Istr, _A0);
		_In(_Istr, _B0);
		_Par._Init(_A0, _B0);
		return (_Istr);
	}

	template<class _Elem,
	class _Traits>
		basic_ostream<_Elem, _Traits>& _Write(
		basic_ostream<_Elem, _Traits>& _Ostr) const
	{	// write state to _Ostr
		_Out(_Ostr, _Par._Ax);
		_Out(_Ostr, _Par._Bx);
		return (_Ostr);
	}

private:
	template<class _Engine>
	result_type _Eval(_Engine& _Eng, const param_type& _Par0) const
	{	// generate pseudo-random value
		_Ty _Px = _NRAND(_Eng, _Ty);
		return (_Par0._Ax + _Par0._Bx*_CSTD log(-_CSTD log(1 - _Px)));
		//return (_Par0._Ax - _Par0._Bx * _CSTD log(-_CSTD log(_Px)));
	}

	param_type _Par;
};

_STD_END
#pragma pop_macro("new")
#pragma warning(pop)
#pragma pack(pop)
