// --------------------------------------------------------------------------
//    This part of code borrows from OpenMS(Open-Source Mass Spectrometry)
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
#include<sstream>
#include<vector>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax2/Attributes.hpp>
using namespace std;
using namespace xercesc;

/**
@brief Base class for XML handlers.
*/ 

class  XMLHandler :
	public xercesc::DefaultHandler
{
public:
	//Action to set the current mode (for error messages)
	enum ActionMode
	{
		LOAD,               ///< Loading a file
		STORE               ///< Storing a file
	};

	// Default constructor
	XMLHandler(const string & filename);
	/// Destructor
	virtual ~XMLHandler();

	// Release internal memory used for parsing (call
	void reset();

	/**
	@name Reimplemented XERCES-C error handlers
	These methods forward the error message to our own error handlers below.
	*/

	// Parsing method for character data
	virtual void characters(const XMLCh * const chars, const XMLSize_t length);
	// Parsing method for opening tags
	virtual void startElement(const XMLCh * const uri, const XMLCh * const localname, const XMLCh * const qname, const xercesc::Attributes & attrs);
	// Parsing method for closing tags
	virtual void endElement(const XMLCh * const uri, const XMLCh * const localname, const XMLCh * const qname);


protected:
	// Error message of the last error
	mutable string error_message_;

	mutable std::vector<XMLCh *> xml_strings_;
	mutable std::vector<char *> c_strings_;
	// File name
	string file_;

	// Schema version
	//string version_;

	/**
	@brief Stack of open XML tags

	This member is used only in those XML parsers that need this information.
	*/

	std::vector<string> open_tags_;

	/// Returns if two xerces strings are equal
	inline bool equal_(const XMLCh * a, const XMLCh * b)
	{
		return xercesc::XMLString::compareString(a, b) == 0;
	}

	//@name General MetaInfo handling (for idXML, featureXML, consensusXML)
	//@{

	//Writes the content of MetaInfoInterface to the file
	//void writeUserParam_(const string & tag_name, std::ostream & os, const MetaInfoInterface & meta, UInt indent) const;

	//@}

	///@name controlled vocabulary handling methods
	//@{

	// Array of CV term lists (one sublist denotes one term and it's children)
	std::vector<std::vector<string> > cv_terms_;

	// Converts @p term to the index of the term in the cv_terms_ entry @p section
	// If the term is not found, @p result_on_error is returned (0 by default)
	inline int cvStringToEnum_(const size_t section, const string & term, const char * message, const int result_on_error = 0)
	{
		//OPENMS_PRECONDITION(section < cv_terms_.size(), "cvStringToEnum_: Index overflow (section number too large)");

		std::vector<string>::const_iterator it = std::find(cv_terms_[section].begin(), cv_terms_[section].end(), term);
		if (it != cv_terms_[section].end())
		{
			return it - cv_terms_[section].begin();
		}
		else
		{
			cout << "Unexpected CV entry " << term << endl;
			return result_on_error;
		}
	}

	//@}
	///@name String conversion
	//@{
	/// Conversion of a string to an integer value
	inline int asInt_(const string & in)
	{
		int res = 0;
		std::stringstream ss(in.c_str());
		if (!(ss >> res))
		{
			cout << "Error:\t When convert string to int\n";
		}

		return res;
	}

	/// Conversion of a Xerces string to an integer value
	inline int asInt_(const XMLCh * in)
	{
		return xercesc::XMLString::parseInt(in);
	}

	/// Conversion of a string to an unsigned integer value
	inline  unsigned int  asUInt_(const string & in)
	{
		unsigned int  res = 0;
		int tmp;
		std::stringstream ss(in.c_str());
		if (!(ss >> tmp))
		{
			cout << "Error:\t When convert string to int\n";
		}
		res = (unsigned int)tmp;
		return res;
	}

	/**
	@brief Conversion of a string to a boolean value
	'true', 'false', '1' and '0' are accpeted.
	@n For all other values a parse error is produced.
	*/

	inline bool asBool_(const string & in)
	{
		if (in == "true" || in == "TRUE" || in == "True" || in == "1")
		{
			return true;
		}
		else if (in == "false" || in == "FALSE" || in == "False" || in == "0")
		{
			return false;
		}
		else
		{
			cout << "Boolean conversion error " << endl;
		}
		return false;
	}

	///@name Accessing attributes
	//@{
	// Converts an attribute to a String
	inline char * attributeAsString_(const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0) cout << "Required attribute '" << name << "' not present!\n";//fatalError(LOAD, String("Required attribute '") + name + "' not present!");
		char * pChar = xercesc::XMLString::transcode(val);
		c_strings_.push_back(pChar);
		return pChar;
	}

	// Converts an attribute to a Int
	inline int attributeAsInt_(const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0) cout << "Required attribute '" << name << "' not present!\n";//fatalError(LOAD, String("Required attribute '") + name + "' not present!");

		const char *pchar = "charge";
		if (!strcmp(name, pchar))
		{
			char * cTemp = xercesc::XMLString::transcode(val);
			if (cTemp[1] != '\0')
				cTemp[1] = '\0';
			return atoi(cTemp);
		}
		else
			return xercesc::XMLString::parseInt(val);

	}

	// Converts an attribute to a DoubleReal
	inline double attributeAsDouble_(const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0) cout << "Required attribute '" << name << "' not present!\n";//fatalError(LOAD, String("Required attribute '") + name + "' not present!");
		char * pChar = xercesc::XMLString::transcode(val);
		c_strings_.push_back(pChar);

		return atof(pChar);
	}

	/**
	@brief Assigns the attribute content to the String @a value if the attribute is present
	@return if the attribute was present
	*/

	inline bool optionalAttributeAsString_(string & value, const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0)
		{
			cout << "Required attribute '" << name << "' not present!\n";//fatalError(LOAD, String("Required attribute '") + name + "' not present!");
			return false;
		}
		char * pChar = xercesc::XMLString::transcode(val);
		c_strings_.push_back(pChar);
		value = pChar;
		return true;
	}

	/**
	@brief Assigns the attribute content to the Int @a value if the attribute is present
	@return if the attribute was present
	*/
	inline bool optionalAttributeAsInt_(int & value, const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0)
		{
			cout << "Required attribute '" << name << "' not present!\n";//fatalError(LOAD, String("Required attribute '") + name + "' not present!");
			return false;
		}

		const char *pchar = "charge";
		if (!strcmp(name, pchar))
		{
			char * cTemp = xercesc::XMLString::transcode(val);
			if (cTemp[1] != '\0')
				cTemp[1] = '\0';
			value = atoi(cTemp);
		}
		else
			value = xercesc::XMLString::parseInt(val);
		return true;
	}

	/**
	@brief Assigns the attribute content to the UInt @a value if the attribute is present
	@return if the attribute was present
	*/

	inline bool optionalAttributeAsUInt_(unsigned int & value, const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0)
		{
			cout << "Required attribute '" << name << "' not present!\n";//fatalError(LOAD, String("Required attribute '") + name + "' not present!");
			return false;
		}

		const char *pchar = "charge";
		if (!strcmp(name, pchar))
		{
			char * cTemp = xercesc::XMLString::transcode(val);
			if (cTemp[1] != '\0')
				cTemp[1] = '\0';
			value = atoi(cTemp);
		}
		else
			value = xercesc::XMLString::parseInt(val);
		return true;
	}

	/**
	@brief Assigns the attribute content to the DoubleReal @a value if the attribute is present
	@return if the attribute was present
	*/

	inline bool optionalAttributeAsDouble_(double & value, const xercesc::Attributes & a, const char * name) const
	{
		XMLCh * result = xercesc::XMLString::transcode(name);
		xml_strings_.push_back(result);
		const XMLCh * val = a.getValue(result);
		if (val == 0)
		{
			cout << "Required attribute '" << name << "' not present!\n";
			return false;
		}
		char * pChar = xercesc::XMLString::transcode(val);
		c_strings_.push_back(pChar);

		value = atof(pChar);
		return true;
	}

private:
	/// Not implemented
	XMLHandler();

	inline string expectList_(const char * str) const
	{
		string tmp(str);
		bool bhasPrefix = false;
		bool bhasSuffix = false;
		string strPrefix = "[";
		string strSufix = "]";

		if (tmp.compare(0, strPrefix.size(), strPrefix) == 0)
			bhasPrefix = true;
		if (tmp.compare(tmp.size() - strSufix.size(), strSufix.size(), strSufix) == 0)
			bhasSuffix = true;
		if (!(bhasPrefix&&bhasSuffix))
		{
			cout << "List argument is not a string representation of a list!\n";
		}
		return tmp;
	}

};

#endif // OPENMS_FORMAT_HANDLERS_XMLHANDLER_H
