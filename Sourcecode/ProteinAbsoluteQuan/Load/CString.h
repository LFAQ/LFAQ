// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_STRING_H
#define OPENMS_DATASTRUCTURES_STRING_H

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
using namespace std;

/**
@brief A more convenient string class.

It is based on std::string but adds a lot of methods for convenience.

@ingroup Datastructures
*/
class String :
	public std::string
{
public:

	/// Empty string for comparisons
	static const String EMPTY;

	/** @name Type definitions
	*/
	//@{
	/// Iterator
	typedef iterator    Iterator;
	/// Const Iterator
	typedef const_iterator  ConstIterator;
	/// Reverse Iterator
	typedef reverse_iterator    ReverseIterator;
	/// Const reverse Iterator
	typedef const_reverse_iterator  ConstReverseIterator;
	/// UInt type
	typedef size_type   SizeType;

	/// How to handle embedded quotes when quoting strings
	enum QuotingMethod { NONE, ESCAPE, DOUBLE };

	//@}

	/**	@name Constructors
	*/
	//@{
	/// Default constructor
	String();
	/// Constructor from std::string
	String(const std::string & s);
	String(const char * s);
	/// Constructor from a char
	String(const char c);
	/// Constructor from char* (only @p length chracters)
	String(const char * s, SizeType length);
	/// Constructor from char (repeats the char @p len times)
	String(size_t len, char c);
	/// Constructor from a char range
	template <class InputIterator>
	String(InputIterator first, InputIterator last) :
		std::string(first, last)
	{

	}

	/// Constructor from an integer
	String(int i);
	/// Constructor from an unsigned integer
	String(unsigned int i);
	/// Constructor from an integer
	String(short int i);
	/// Constructor from an unsigned integer
	String(short unsigned int i);
	/// Constructor from an integer
	String(long int i);
	/// Constructor from an unsigned integer
	String(long unsigned int i);
	/// Constructor from an unsigned integer
	String(long long unsigned int i);
	/// Constructor from an unsigned integer
	String(long long signed int i);
	/// Constructor from float
	String(float f);
	/// Constructor from double
	String(double d);
	/// Constructor from long double
	String(long double ld);
	/// Constructor from DataValue (casted to String)

	//@}

	/** @name Predicates
	*/
	//@{
	/// true if String begins with @p string, false otherwise
	bool hasPrefix(const String & string) const;

	/// true if String ends with @p string, false otherwise
	bool hasSuffix(const String & string) const;

	/// true if String contains the @p string, false otherwise
	bool hasSubstring(const String & string) const;

	/// true if String contains the @p byte, false otherwise
	bool has(unsigned char byte) const;
	//@}


	/** @name Accessors
	*/
	//@{
	/**
	@brief returns the prefix of length @p length

	@exception Exception::IndexOverflow is thrown if @p length is bigger than the size
	*/
	String prefix(SizeType length) const;

	/**
	@brief returns the suffix of length @p length

	@exception Exception::IndexOverflow is thrown if @p length is bigger than the size
	*/
	String suffix(SizeType length) const;

	/**
	@brief returns the prefix of length @p length

	@exception Exception::IndexUnderflow is thrown if @p length is smaller than zero
	@exception Exception::IndexOverflow is thrown if @p length is bigger than the size
	*/
	String prefix(int length) const;

	/**
	@brief returns the suffix of length @p length

	@exception Exception::IndexUnderflow is thrown if @p length is smaller than zero
	@exception Exception::IndexOverflowis thrown if @p length is bigger than the size
	*/
	String suffix(int length) const;

	/**
	@brief returns the prefix up to the first occurence of char @p delim (excluding it)

	@exception Exception::ElementNotFound is thrown if @p delim is not found
	*/
	String prefix(char delim) const;

	/**
	@brief returns the suffix up to the last occurence of char @p delim (excluding it)

	@exception Exception::ElementNotFound is thrown if @p delim is not found
	*/
	String suffix(char delim) const;

	/**
	@brief Wrapper for the STL substr() method. Returns a String object with its contents initialized to a substring of the current object.

	@param pos Position of a character in the current string object to be used as starting character for the substring.
	If the @p pos is past the end of the string, it is set to the end of the string.

	@param n Length of the substring.
	If this value would make the substring to span past the end of the current string content, only those characters until the end of the string are used.
	npos is a static member constant value with the greatest possible value for an element of type size_t, therefore, when this value is used, all the
	characters between pos and the end of the string are used as the initialization substring.

	*/
	String substr(size_t pos = 0, size_t n = npos) const;

	/**
	@brief Returns a substring where @p n characters were removed from the end of the string.

	If @p n is greater than size(), the result is an empty string.

	@param n Number of characters that will be removed from the end of the string.
	*/
	String chop(size_t n) const;

	//@}


	/**
	@name Mutators

	All these methods return a reference to the string in order to make them chainable
	*/
	//@{
	/// inverts the direction of the string
	String & reverse();

	/// removes whitespaces (space, tab, line feed, carriage return) at the beginning and the end of the string
	String & trim();

	/**
	@brief Wraps the string in quotation marks

	The quotation mark can be specified by parameter @p q (typically single or double quote); embedded quotation marks are handled according to @p method by backslash-escaping, doubling, or not at all.

	@see unquote()
	*/
	String & quote(char q = '"', QuotingMethod method = ESCAPE);

	/**
	@brief Reverses changes made by the @p quote method

	Removes surrounding quotation marks (given by parameter @p q); handles embedded quotes according to @p method.

	@exception Exception::ConversionError is thrown if the string does not have the format produced by @p quote

	@see quote()
	*/
	String & unquote(char q = '"', QuotingMethod method = ESCAPE);

	/// merges subsequent whitespaces to one blank character
	String & simplify();

	///Adds @p c on the left side until the size of the string is @p size
	String & fillLeft(char c, int size);

	///Adds @p c on the right side until the size of the string is @p size
	String & fillRight(char c, int size);

	///Converts the string to uppercase
	String & toUpper();

	///Converts the string to lowercase
	String & toLower();

	///Converts the first letter of the string to uppercase
	String & firstToUpper();

	///Replaces all occurences of the character @p from by the character @p to.
	String & substitute(char from, char to);

	///Replaces all occurences of the string @p from by the string @p to.
	String & substitute(const String & from, const String & to);

	///Remove all occurences of the character @p what.
	String & remove(char what);

	///Makes sure the string ends with the character @p end
	String & ensureLastChar(char end);

	///removes whitespaces (space, tab, line feed, carriage return)
	String & removeWhitespaces();
	//@}

	/** @name Converters
	*/
	//@{

	/**
	@brief Conversion to int

	This method extracts only the integral part of the string.
	If you want the result rounded, use toFloat() and round the result.

	@exception Exception::ConversionError is thrown if the string could not be converted to int
	*/
	int toInt() const;

	/**
	@brief Conversion to float

	@exception Exception::ConversionError is thrown if the string could not be converted to float
	*/
	float toFloat() const;

	/**
	@brief Conversion to double

	@exception Exception::ConversionError is thrown if the string could not be converted to double
	*/
	double toDouble() const;

	/// Conversion to Qt QString
	//	 QString toQString() const;

	//@}

	/** @name Sum operator overloads
	*/
	//@{
	/// Sum operator for an integer
	String operator+(int i) const;
	/// Sum operator for an unsigned integer
	String operator+(unsigned int i) const;
	/// Sum operator for an integer
	String operator+(short int i) const;
	/// Sum operator for an unsigned integer
	String operator+(short unsigned int i) const;
	/// Sum operator for an integer
	String operator+(long int i) const;
	/// Sum operator for an unsigned integer
	String operator+(long unsigned int i) const;
	/// Sum operator for an unsigned integer
	String operator+(long long unsigned int i) const;
	/// Sum operator for float
	String operator+(float f) const;
	/// Sum operator for double
	String operator+(double d) const;
	/// Sum operator for long double
	String operator+(long double ld) const;
	/// Sum operator for char
	String operator+(char c) const;
	/// Sum operator for char*
	String operator+(const char * s) const;
	/// Sum operatr for String
	String operator+(const String & s) const;
	/// Sum operator for std::string
	String operator+(const std::string & s) const;
	//@}

	/** @name Append operator overloads
	*/
	//@{
	/// Sum operator for an integer
	String & operator+=(int i);
	/// Sum operator for an unsigned integer
	String & operator+=(unsigned int i);
	/// Sum operator for an integer
	String & operator+=(short int i);
	/// Sum operator for an unsigned integer
	String & operator+=(short unsigned int i);
	/// Sum operator for an integer
	String & operator+=(long int i);
	/// Sum operator for an unsigned integer
	String & operator+=(long unsigned int i);
	/// Sum operator for an unsigned integer
	String & operator+=(long long unsigned int i);
	/// Sum operator for float
	String & operator+=(float f);
	/// Sum operator for double
	String & operator+=(double d);
	/// Sum operator for long double
	String & operator+=(long double d);
	/// Sum operator for char
	String & operator+=(char c);
	/// Sum operator for char*
	String & operator+=(const char * s);
	/// Sum operatr for String
	String & operator+=(const String & s);
	/// Sum operator for std::string
	String & operator+=(const std::string & s);
	//@}

	///returns a random string of the given length. It consists of [0-9a-zA-Z]
	static String random(int length);

	///returns a string for @p d with exactly @p n decimal places
	// static String number(double d, int n);
	/**
	@brief Returns a string with at maximum @p n characters for @p d

	If @p d is larger, scientific notation is used.
	*/
	static String numberLength(double d, int n);


	/**
	@brief Splits a string into @p substrings using @p splitter as delimiter

	If @p splitter is not found, the whole string is put into @p substrings.
	If @p splitter is empty, the string is split into individual characters.
	If the invoking string is empty, @p substrings will also be empty.

	@p quote_protect (default: false) can be used to split only between quoted
	blocks e.g. ' "a string" , "another string with , in it" '
	results in only two substrings (with double quotation marks @em removed).
	Every returned substring is trimmed and then (if present) has surrounding quotation marks removed.

	@return @e true if one or more splits occurred, @e false otherwise

	@see concatenate().
	*/
	bool split(const char splitter, std::vector<String> & substrings, bool quote_protect = false) const;

	/**
	@brief Splits a string into @p substrings using @p splitter (the whole string) as delimiter

	If @p splitter is not found,  the whole string is put into @p substrings.
	If @p splitter is empty, the string is split into individual characters.
	If the invoking string is empty, @p substrings will also be empty.

	@return @e true if one or more splits occurred, @e false otherwise

	@see concatenate().
	*/
	bool split(const String & splitter, std::vector<String> & substrings) const;

	/**
	@brief Splits a string into @p substrings using @p splitter (the whole string) as delimiter, but does not split within quoted substrings

	A "quoted substring" has the format as produced by @p quote(q, method), where @p q is the quoting character and @p method defines the handling of embedded quotes. Substrings will not be "unquoted" or otherwise processed.

	If @p splitter is not found,  the whole string is put into @p substrings.
	If @p splitter or the invoking string is empty, @p substrings will also be empty.

	@return @e true if one or more splits occurred, @e false otherwise

	@exception Exception::ConversionError is thrown if quotation marks are not balanced

	@see concatenate(), quote().
	*/
	bool split_quoted(const String & splitter, std::vector<String> & substrings,
		char q = '"', QuotingMethod method = ESCAPE) const;

	/**
	@brief Concatenates all elements from @p first to @p last-1 and inserts @p glue between the elements

	@see split().
	*/
	template <class StringIterator>
	void concatenate(StringIterator first, StringIterator last, const String & glue = "")
	{
		//empty container
		if (first == last)
		{
			std::string::clear();
			return;
		}

		std::string::operator=(*first);
		for (StringIterator it = ++first; it != last; ++it)
		{
			std::string::operator+=(glue + (*it));
		}
	}

};



#endif // OPENMS_DATASTRUCTURES_STRING_H
