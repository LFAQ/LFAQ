
#include <set>
#include <utility>
#include"BasicClass.h"
using namespace std;

/**
@brief Representation of a controlled vocabulary.

This representation only contains the information used for parsing and validation.
All other lines are stored in the @em unparsed member of the the CVTerm struct.

@ingroup Format
*/
class  ControlledVocabulary
{
	friend  std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv);

public:
	/// Representation of a CV term
	struct  CVTerm
	{
		/// define xsd types allowed in cv term to specify their value-type
		enum XRefType
		{
			XSD_STRING = 0, // xsd:string A string
			XSD_INTEGER, // xsd:integer Any integer
			XSD_DECIMAL, // xsd:decimal Any real number
			XSD_NEGATIVE_INTEGER, // xsd:negativeInteger Any negative integer
			XSD_POSITIVE_INTEGER, // xsd:positiveInteger Any integer > 0
			XSD_NON_NEGATIVE_INTEGER, // xsd:nonNegativeInteger Any integer >= 0
			XSD_NON_POSITIVE_INTEGER, // xsd:nonPositiveInteger Any integer < 0
			XSD_BOOLEAN, // xsd:boolean True or false
			XSD_DATE, // xsd:date An XML-Schema date
			XSD_ANYURI, // xsd:anyURI uniform resource identifier
			NONE
		};

		static string getXRefTypeName(XRefType type);

		string name; ///< Text name
		string id; ///< Identifier
		std::set<string> parents; ///< The parent IDs
		std::set<string> children; ///< The child IDs
		bool obsolete; ///< Flag that indicates of the term is obsolete
		string description; ///< Term description
		vector<string>synonyms; ///< List of synonyms
		vector<string>unparsed; ///< Unparsed lines from the definition file
		XRefType xref_type; ///< xref value-type for the CV-term
		vector<string>xref_binary; ///< xref binary-data-type for the CV-term (list of all allowed data value types for the current binary data array)
		std::set<string> units; ///< unit accession ids, defined by relationship has units

		///Default constructor
		CVTerm();

		CVTerm(const CVTerm& rhs);

		CVTerm& operator=(const CVTerm& rhs);

		/// get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
		string toXMLString(const string& ref, const string& value = string("")) const;

	};

	/// Constructor
	ControlledVocabulary();

	///Destructor
	virtual ~ControlledVocabulary();

	/// Returns the CV name (set in the load method)
	const string& name() const;

	/**
	@brief Loads the CV from an OBO file

	@exception Exception::FileNotFound is thrown if the file could not be opened
	@exception Exception::ParseError is thrown if an error occurs during parsing
	*/
	void loadFromOBO(const string& name, const string& filename);

	/// Returns true if the term is in the CV. Returns false otherwise.
	bool exists(const string& id) const;

	/// Returns true if a term with the given name is in the CV. Returns false otherwise.
	bool hasTermWithName(const string& name) const;

	/**
	@brief Returns a term specified by ID
	@exception Exception::InvalidValue is thrown if the term is not present
	*/
	const CVTerm& getTerm(const string& id) const;

	/**
	@brief Returns a term specified by name

	@exception Exception::InvalidValue is thrown if the term is not present
	*/
	const CVTerm& getTermByName(const string& name, const string& desc = "") const;

	/// returns all the terms stored in the CV
	const map<string, CVTerm>& getTerms() const;

	/**
	@brief Writes all child terms recursively into terms
	If parent has child this method writes them recursively into the term object
	@exception Exception::InvalidValue is thrown if the term is not present
	*/
	void getAllChildTerms(std::set<string>& terms, const string& parent) const;

	/**
	@brief Returns if @p child is a child of @p parent
	@exception Exception::InvalidValue is thrown if one of the terms is not present
	*/
	bool isChildOf(const string& child, const string& parent) const;

protected:
	/**
	@brief checks if a name corresponds to an id
	If the term is not known, 'true' is returned!
	*/
	bool checkName_(const string& id, const string& name, bool ignore_case = true);

	///Map from ID to CVTerm
	map<string, CVTerm> terms_;
	///Map from name to id
	map<string, string> namesToIds_;
	///Name set in the load method
	string name_;
};

///Print the contents to a stream.
std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv);


