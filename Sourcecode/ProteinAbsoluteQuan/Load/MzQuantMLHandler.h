#include "XMLHandler.h"
#include "ControlledVocabulary.h"
#include "QuantInformation.h"
#include  <xercesc/util/XMLString.hpp>  

using namespace std;
/**
@brief XML handler for mzQuantMLFile
@note Do not use this class. It is only needed in CloadmzQuantML.
*/

class  mzQuantMLHandler :
	public XMLHandler
{
public:
	/**@name Constructors and destructor */
	mzQuantMLHandler(CQuantInformation& qi, const string & filename);

	/// Destructor
	virtual ~mzQuantMLHandler();
	//@}

	// Docu in base class
	virtual void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname);

	// Docu in base class
	virtual void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes);

	// Docu in base class
	virtual void characters(const XMLCh * const chars, const XMLSize_t length);

protected:
	/// Controlled vocabulary (hopefully the psi-pi from OpenMS/share/OpenMS/CV/psi-pi.obo)
	ControlledVocabulary cv_;
	string tag_;
	CQuantInformation *qi_;

private:
	mzQuantMLHandler();
	mzQuantMLHandler(const mzQuantMLHandler & rhs);
	mzQuantMLHandler & operator=(const mzQuantMLHandler & rhs);
	string current_id_;
	std::vector<string> current_col_types_;
	std::vector<double> current_row_;
	std::set<string> to_ignore;
};
