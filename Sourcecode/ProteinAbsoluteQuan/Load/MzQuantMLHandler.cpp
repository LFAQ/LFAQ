#include "stdafx.h"
#include "mzQuantMLHandler.h"


mzQuantMLHandler::mzQuantMLHandler(CQuantInformation &qi, const string & filename) :
XMLHandler(filename), qi_(&qi)
{
	cv_.loadFromOBO("MS", "./CV/psi-ms.obo"); //TODO unimod -> then automatise CVList writing
}

mzQuantMLHandler::~mzQuantMLHandler()
{
}

void mzQuantMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
{
	char * result = xercesc::XMLString::transcode(qname);
	c_strings_.push_back(result);
	tag_ = result;
	open_tags_.push_back(tag_);
	static set<string> to_ignore;

	if (to_ignore.empty())
	{
		to_ignore.insert("CvList"); // for now static set of obos.
		to_ignore.insert("Cv"); // for now static set of obos.
		to_ignore.insert("ProteinGroupList"); // for now no proteins or groups
		to_ignore.insert("StudyVariableList"); // We can't deal with these right now, but that is coming
		to_ignore.insert("StudyVariable"); 
		to_ignore.insert("Assay_refs"); 
		to_ignore.insert("FeatureList"); // we only need to see the features and datamatrices rows
		to_ignore.insert("AssayList"); // we only need to see the assays
		to_ignore.insert("DataProcessingList"); // we only need to see the DataProcessings
		to_ignore.insert("SoftwareList"); // we only need to see the Softwares
		to_ignore.insert("InputFiles"); // we only need to see the Files
		to_ignore.insert("Label"); // we only need to see the Modifications
		to_ignore.insert("DataType"); // we only need to see the Modifications
		to_ignore.insert("DataMatrix"); // we only need to see the inside characters
	}

	if (to_ignore.find(tag_) != to_ignore.end())
	{
		return;
	}

	//determine parent tag
	string parent_tag;
	if (open_tags_.size() > 1)
	{
		parent_tag = *(open_tags_.end() - 2);
	}
	string parent_parent_tag;
	if (open_tags_.size() > 2)
	{
		parent_parent_tag = *(open_tags_.end() - 3);
	}

	static const XMLCh* s_value = xercesc::XMLString::transcode("value");
	static const XMLCh* s_type = xercesc::XMLString::transcode("type");
	static const XMLCh* s_name = xercesc::XMLString::transcode("name");
	static const XMLCh* s_unit_accession = xercesc::XMLString::transcode("unitAccession");
	static const XMLCh* s_cv_ref = xercesc::XMLString::transcode("cvRef");
	static const XMLCh* s_accession = xercesc::XMLString::transcode("accession");

	if (tag_ == "ColumnIndex")
	{
		qi_->vecExperimentNames.clear();
	}
	else if (tag_ == "AssayQuantLayer")
	{
		current_col_types_.clear();
	}
	else if (tag_ == "PeptideConsensusList")
	{
	}
	else if (tag_ == "Protein")
	{
		if (parent_tag == "ProteinList")
		{
			qi_->strProteinAccession = attributeAsString_(attributes, "accession");
		}

	}
	else if (tag_ == "mzQuantML")
	{
		// handle version and experiment type
	}
	else if (tag_ == "AnalysisSummary")
	{
		// handle version and experiment type
	}
	else if (tag_ == "Row")
	{
		qi_->StrPeptideId = attributeAsString_(attributes, "object_ref");
		current_row_.clear();
	}
	else if (tag_ == "PeptideConsensus")
	{
		qi_->StrPeptideId = attributeAsString_(attributes, "id");
	}
}

void mzQuantMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
{
	//if there is data between the element tags - !attention if element is derived from a xsd:list type, each list entry is a charecters call
	if (tag_ == "PeptideSequence")
	{
		char * result = xercesc::XMLString::transcode(chars);
		c_strings_.push_back(result);
		current_id_ = result;
	}
	else if (tag_ == "Row")
	{
		char * result = xercesc::XMLString::transcode(chars);
		c_strings_.push_back(result);

		string r = result;
		r = stringTrim(r);
		if (!r.empty()) // always two notifications for a row, only the first one contains chars - dunno why
		{
			std::vector<string> splits;
			stringSplit(r, " ", splits);
			for (std::vector<string>::iterator it = splits.begin(); it != splits.end(); ++it)
			{
				current_row_.push_back(stod(*it));
			}
		}
	}
	else if (tag_ == "PeptideConsensus_refs")
	{
		char * result = xercesc::XMLString::transcode(chars);
		c_strings_.push_back(result);

		string r = result;
		r = stringTrim(r);

		if (!r.empty()) // always two notifications for a row, only the first one contains chars - dunno why
		{
			std::vector<string> splits;
			stringSplit(r, " ", splits);
			for (std::vector<string>::iterator it = splits.begin(); it != splits.end(); ++it)
			{
				qi_->vecPeptideIds.push_back(*it);
			}
		}
	}
	else if (tag_ == "ColumnIndex")
	{
		//overwrites current_col_types_ with the ratio_refs or the assay_refs
		char * result = xercesc::XMLString::transcode(chars);
		c_strings_.push_back(result);

		string r = result;
		//clear must have happened earlyer in QuantLayer tag
		r = stringTrim(r);
		if (!r.empty()) // always two notifications for a row, only the first one contains chars - dunno why
		{
			stringSplit(r, " ", current_col_types_);
			stringSplit(r, " ", qi_->vecExperimentNames);
		}
	}
}

void mzQuantMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
{
	static set<string> to_ignore;
	if (to_ignore.empty())
	{
		to_ignore.insert("Cv");
	}

	char * result = xercesc::XMLString::transcode(qname);
	c_strings_.push_back(result);

	tag_ = result;
	//determine parent tag
	string parent_tag;
	if (open_tags_.size() > 1)
	{
		parent_tag = *(open_tags_.end() - 2);
	}
	string parent_parent_tag;
	if (open_tags_.size() > 2)
	{
		parent_parent_tag = *(open_tags_.end() - 3);
	}

	string parent_parent_parent_tag;
	if (open_tags_.size() > 3)
	{
		parent_parent_parent_tag = *(open_tags_.end() - 4);
	}
	//close current tag
	open_tags_.pop_back();

	if (to_ignore.find(tag_) != to_ignore.end())
	{
		return;
	}

	// no ProcessingMethod endElement action so each userParam under Dataprocessing will be one processingaction 
	// no other way for core-lib compability yet
	if (tag_ == "DataProcessing")
	{
		return;
	}
	else if (tag_ == "Protein")
	{
		qi_->mapProteinAndPeptideIDs.insert(pair<string, vector<string>>(qi_->strProteinAccession, qi_->vecPeptideIds));
		qi_->vecPeptideIds.clear();
	}
	else if (tag_ == "DataProcessingList")
	{
	}
	else if (tag_ == "ColumnDefinition")
	{
		//TODO check all current_col_types_[] are not empty
	}
	else if (tag_ == "Row")
	{
		if (parent_parent_tag == "AssayQuantLayer" &&current_col_types_.size() != current_row_.size())
		{
			cout << "Unknown / unmatching row content in Row element of " << parent_tag << "'.";
			return;
		}

		if (parent_parent_tag == "AssayQuantLayer"&&parent_parent_parent_tag == "PeptideConsensusList")
		{
			map<string, double> mapExpeimentAndIntensity;

			for (int i = 0; i < current_col_types_.size(); i++)
			{
				mapExpeimentAndIntensity.insert(pair<string, double>(current_col_types_.at(i), current_row_.at(i)));
			}
			qi_->mapPeptideAndExperimenIntensity.insert(pair<string, map<string, double>>(qi_->StrPeptideId, mapExpeimentAndIntensity));
		}
	}
	else if (tag_ == "PeptideSequence")
	{
		qi_->mapPeptideIdAndSequence.insert(pair<string, string>(qi_->StrPeptideId, current_id_));
	}
}
