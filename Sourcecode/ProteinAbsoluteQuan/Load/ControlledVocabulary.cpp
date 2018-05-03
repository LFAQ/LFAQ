#include "stdafx.h"
#include "ControlledVocabulary.h"


ControlledVocabulary::CVTerm::CVTerm() :
name(),
id(),
parents(),
children(),
obsolete(false),
description(),
synonyms(),
unparsed(),
xref_type(NONE),
xref_binary()
{
}

ControlledVocabulary::CVTerm::CVTerm(const CVTerm& rhs) :
name(rhs.name),
id(rhs.id),
parents(rhs.parents),
children(rhs.children),
obsolete(rhs.obsolete),
description(rhs.description),
synonyms(rhs.synonyms),
unparsed(rhs.unparsed),
xref_type(rhs.xref_type),
xref_binary(rhs.xref_binary),
units(rhs.units)
{
}

ControlledVocabulary::CVTerm& ControlledVocabulary::CVTerm::operator=(const CVTerm& rhs)
{
	if (this != &rhs)
	{
		name = rhs.name;
		id = rhs.id;
		parents = rhs.parents;
		children = rhs.children;
		obsolete = rhs.obsolete;
		description = rhs.description;
		synonyms = rhs.synonyms;
		unparsed = rhs.unparsed;
		xref_type = rhs.xref_type;
		xref_binary = rhs.xref_binary;
		units = rhs.units;
	}
	return *this;
}

string ControlledVocabulary::CVTerm::getXRefTypeName(XRefType type)
{
	switch (type)
	{
	case XSD_STRING: return "xsd:string";

	case XSD_INTEGER: return "xsd:integer";

	case XSD_DECIMAL: return "xsd:decimal";

	case XSD_NEGATIVE_INTEGER: return "xsd:negativeInteger";

	case XSD_POSITIVE_INTEGER: return "xsd:positiveInteger";

	case XSD_NON_NEGATIVE_INTEGER: return "xsd:nonNegativeInteger";

	case XSD_NON_POSITIVE_INTEGER: return "xsd:nonPositiveInteger";

	case XSD_BOOLEAN: return "xsd:boolean";

	case XSD_DATE: return "xsd:date";

	case XSD_ANYURI: return "xsd:anyURI";

	default: return "none";
	}
	return "";
}

string ControlledVocabulary::CVTerm::toXMLString(const string& ref, const string& value) const
{
	string s = "<cvParam accession=\"" + id + "\" cvRef=\"" + ref + "\" name=\"" + name;
	if (!value.empty())
	{
		s += "\" value=\"" + value;
	}
	s += "\"/>";
	return s;
	//~ TODO: handle unknown cvparams in ControlledVocabulary to get same formatting but more userdefined interface
}

ControlledVocabulary::ControlledVocabulary() :
terms_(),
name_("")
{

}

ControlledVocabulary::~ControlledVocabulary()
{

}

void ControlledVocabulary::loadFromOBO(const string& name, const string& filename)
{
	bool in_term = false;
	name_ = name;

	ifstream is(filename.c_str());
	if (!is)
	{
		cout << "Cannot find file " << filename << endl;
        flog.mf_Input("Cannot find file " + filename + "\n");
        flog.mf_Destroy();
        exit(-1);
	}

	string line, line_wo_spaces;
	CVTerm term;

	//parse file
	while (getline(is, line, '\n'))
	{
		line = stringTrim(line);
		line_wo_spaces = line;
		line_wo_spaces = removestringWhitespaces(line_wo_spaces);

		//do nothing for empty lines
		if (line == "")
			continue;

		//********************************************************************************
		//stanza line
		if (line_wo_spaces[0] == '[')
		{
			if (stringtoLower(line_wo_spaces) == "[term]") //new term
			{
				in_term = true;
				if (term.id != "") //store last term
				{
					terms_[term.id] = term;
				}
				//clear temporary term members
				term = CVTerm();
			}
			// other stanza => not in a term
			else
			{
				in_term = false;
			}
		}
		//********************************************************************************
		//data line
		else if (in_term)
		{
			if (bPrefix(line_wo_spaces,"id:"))
			{
				term.id = stringTrim(line.substr(line.find(':') + 1));
			}
			else if (bPrefix(line_wo_spaces,"name:"))
			{
				term.name = stringTrim(line.substr(line.find(':') + 1));
			}
			else if (bPrefix(line_wo_spaces,"is_a:"))
			{
				/*if (line.has('!'))*/
				if (line.find('!')!=string::npos)
				{
					string parent_id =stringTrim(PrefixOfstring(line.substr(line.find(':') + 1),'!'));
					term.parents.insert(parent_id);

					//check if the parent term name is correct
					string parent_name = stringTrim(SuffixOfstring(line,'!'));
					if (!checkName_(parent_id, parent_name))
						cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': parent term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
				}
				else
				{
					term.parents.insert(stringTrim(line.substr(line.find(':') + 1)));
				}
			}
			// brenda tissue special relationships, DRV (derived and part of)
			else if (bPrefix(line_wo_spaces,"relationship:DRV") && name == "brenda")
			{
				if (line.find('!')!=string::npos)
				{
					// e.g. relationship: DRV BTO:0000142 ! brain
					string parent_id = PrefixOfstring(line.substr(line.find("DRV") + 4), ':') + ":" + stringTrim(PrefixOfstring(SuffixOfstring(line,':'), '!'));
					term.parents.insert(parent_id);

					//check if the parent term name is correct
					string parent_name = stringTrim(SuffixOfstring(line,'!'));
					if (!checkName_(parent_id, parent_name))
						cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': DRV relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
				}
				else
				{
					// e.g. relationship: DRV BTO:0000142
					term.parents.insert(PrefixOfstring(line.substr(line.find("DRV") + 4), ':') + ":" + stringTrim(SuffixOfstring(line,':')));
				}
			}
			else if (bPrefix(line_wo_spaces,"relationship:part_of") && name == "brenda")
			{
				if (line.find('!')!=string::npos)
				{
					string parent_id = PrefixOfstring(line.substr(line.find("part_of") + 8), ':') + ":" + stringTrim(PrefixOfstring(SuffixOfstring(line,':'), '!'));
					term.parents.insert(parent_id);

					//check if the parent term name is correct
					string parent_name = stringTrim(SuffixOfstring(line,'!'));
					if (!checkName_(parent_id, parent_name))
					{
						cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': part_of relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << "\n";
						flog.mf_Input("Warning:\t while loading term '" + term.id + "' of CV '" + name_ + "': part_of relationship term name '" + parent_name + "' and id '" + parent_id + "' differ.\n");
					}
				}
				else
				{
					term.parents.insert(PrefixOfstring(line.substr(line.find("part_of") + 8), ':') + ":" + stringTrim(SuffixOfstring(line,':')));
				}
			}
			else if (bPrefix(line_wo_spaces,"relationship:has_units"))
			{
				if (line.find('!')!=string::npos)
				{
					string unit_id = PrefixOfstring(line.substr(line.find("has_units") + 10), ':') + ":" + stringTrim(PrefixOfstring(SuffixOfstring(line,':'), '!'));
					term.units.insert(unit_id);

					//check if the parent term name is correct
					string unit_name = stringTrim(SuffixOfstring(line,'!'));
					if (!checkName_(unit_id, unit_name))
					{
						cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': has_units relationship term name '" \
							<< unit_name << "' and id '" << unit_id << "' differ." << "\n";
						flog.mf_Input("Warning:\t while loading term '" + term.id + "' of CV '" + name_ + "': part_of relationship term name '"\
							+ unit_name + "' and id '" + unit_id + "' differ.\n");
					}
				}
				else
				{
					term.units.insert(PrefixOfstring(line.substr(line.find("has_units") + 10), ':') + ":" + stringTrim(SuffixOfstring(line,':')));
				}
			}
			else if (bPrefix(line_wo_spaces,"def:"))
			{
				string description = line.substr(line.find('"') + 1);
				description = stringTrim(description);
				description = description.substr(0, description.find('"'));
				description = stringTrim(description);
				term.description = description;
			}
			else if (bPrefix(line_wo_spaces,"synonym:"))
			{
				string synonym = line.substr(line.find('"') + 1);
				synonym=stringTrim(synonym);
				synonym = synonym.substr(0, synonym.find('"'));
				synonym=stringTrim(synonym);

				term.synonyms.push_back(synonym);
			}
			else if (line_wo_spaces == "is_obsolete:true")
			{
				term.obsolete = true;
			}
			else if (bPrefix(line_wo_spaces,"xref:value-type") || bPrefix(line_wo_spaces,"xref_analog:value-type"))
			{
				line_wo_spaces.erase(std::remove(line_wo_spaces.begin(), line_wo_spaces.end(), '\\'), line_wo_spaces.end());
				if (line_wo_spaces.find("value-type:xsd:string")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_STRING;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:integer")!=string::npos || line_wo_spaces.find("value-type:xsd:int")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_INTEGER;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:decimal")!=string::npos ||
					line_wo_spaces.find("value-type:xsd:float")!=string::npos ||
					line_wo_spaces.find("value-type:xsd:double")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_DECIMAL;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:negativeInteger")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_NEGATIVE_INTEGER;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:positiveInteger")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_POSITIVE_INTEGER;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:nonNegativeInteger")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_NON_NEGATIVE_INTEGER;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:nonPositiveInteger")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_NON_POSITIVE_INTEGER;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:boolean")!=string::npos || line_wo_spaces.find("value-type:xsd:bool")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_BOOLEAN;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:date")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_DATE;
					continue;
				}
				if (line_wo_spaces.find("value-type:xsd:anyURI")!=string::npos)
				{
					term.xref_type = CVTerm::XSD_ANYURI;
					continue;
				}
				cerr << "ControlledVocabulary: OBOFile: unknown xsd type: " << line_wo_spaces << ", ignoring" << "\n";
			}
			else if (bPrefix(line_wo_spaces,"xref:binary-data-type") || bPrefix(line_wo_spaces,"xref_analog:binary-data-type"))
			{
				line_wo_spaces.erase(std::remove(line_wo_spaces.begin(), line_wo_spaces.end(), '\\'), line_wo_spaces.end());
				//remove description (if present)
				// according to rev1165 of the cv comments are here quoted,
				//see http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo?revision=1.165&view=markup
				if (line_wo_spaces.find('\"')!=string::npos)
				{
					line_wo_spaces = line_wo_spaces.substr(0, line_wo_spaces.find('\"'));
				}
				//trim prefix
				line_wo_spaces = line_wo_spaces.substr(22);
				//trim just to be sure
				line_wo_spaces =stringTrim(line_wo_spaces);
				term.xref_binary.push_back(line_wo_spaces);
			}
			else if (line != "")
			{
				term.unparsed.push_back(line);
			}
		}
	}

	if (term.id != "") //store last term
	{
		terms_[term.id] = term;
	}

	// now build all child terms
	for (map<string, CVTerm>::iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		for (set<string>::const_iterator pit = terms_[it->first].parents.begin(); pit != terms_[it->first].parents.end(); ++pit)
		{
			terms_[*pit].children.insert(it->first);
		}

		map<string, string>::iterator mit = namesToIds_.find(terms_[it->first].name);
		if (mit == namesToIds_.end())
		{
			namesToIds_.insert(pair<string, string>(terms_[it->first].name, it->first));
		}
		else
		{
			//~ TODO that case would be bad do something
			string s = terms_[it->first].name + terms_[it->first].description;
			namesToIds_.insert(pair<string, string>(s, it->first));
		}
	}
}

const ControlledVocabulary::CVTerm& ControlledVocabulary::getTerm(const string& id) const
{
	map<string, CVTerm>::const_iterator it = terms_.find(id);
	if (it == terms_.end())
	{
		cout << "Warning:\tInvalid CV identifier!" << id << endl;
        flog.mf_Input("Warning:\tInvalid CV identifier!" + id + "\n");
	}
	return it->second;
}

const map<string, ControlledVocabulary::CVTerm>& ControlledVocabulary::getTerms() const
{
	return terms_;
}

void ControlledVocabulary::getAllChildTerms(set<string>& terms, const string& parent) const
{
	const set<string>& children = getTerm(parent).children;
	for (set<string>::const_iterator it = children.begin(); it != children.end(); ++it)
	{
		terms.insert(*it);
		//TODO This is not save for cyclic graphs. Are they allowd in CVs?
		getAllChildTerms(terms, *it);
	}
}

const ControlledVocabulary::CVTerm& ControlledVocabulary::getTermByName(const string& name, const string& desc) const
{
	//slow, but Vocabulary is very finite and this method will be called only a few times during write of a ML file using a CV
	map<string, string>::const_iterator it = namesToIds_.find(name);
	if (it == namesToIds_.end())
	{
		if (!desc.empty())
		{
			it = namesToIds_.find(string(name + desc));
			if (it == namesToIds_.end())
			{
				cout << "Invalid CV name!" << name << endl;
			}
		}
		else
		{
			cout << "Invalid CV name!" << name << endl;
		}
	}
	map<string, CVTerm>::const_iterator termsIter;
	termsIter = terms_.find(it->second);

	return termsIter->second;
}

bool ControlledVocabulary::exists(const string& id) const
{
	if (terms_.find(id) != terms_.end())
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool ControlledVocabulary::hasTermWithName(const string& name) const
{
	map<string, string>::const_iterator it = namesToIds_.find(name);
	return it != namesToIds_.end();
}

bool ControlledVocabulary::isChildOf(const string& child, const string& parent) const
{
	const CVTerm& ch = getTerm(child);

	for (set<string>::const_iterator it = ch.parents.begin(); it != ch.parents.end(); ++it)
	{
		//check if it is a direct parent
		if (*it == parent)
		{
			return true;
		}
		//check if it is an indirect parent
		else if (isChildOf(*it, parent))
		{
			return true;
		}
	}

	return false;
}

std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv)
{
	for (map<string, ControlledVocabulary::CVTerm>::const_iterator it = cv.terms_.begin(); it != cv.terms_.end(); ++it)
	{
		os << "[Term]\n";
		os << "id: '" << it->second.id << "'\n";
		os << "name: '" << it->second.name << "'\n";
		for (set<string>::const_iterator it2 = it->second.parents.begin(); it2 != it->second.parents.end(); ++it2)
		{
			cout << "is_a: '" << *it2 << "'\n";
		}
	}
	return os;
}

const string& ControlledVocabulary::name() const
{
	return name_;
}

bool ControlledVocabulary::checkName_(const string& id, const string& name, bool ignore_case)
{
	if (!exists(id))
		return true;

	string parent_name = name;
	string real_parent_name = getTerm(id).name;

	if (ignore_case)
	{
		parent_name=stringtoLower(parent_name);
		real_parent_name =stringtoLower(real_parent_name);
	}

	return real_parent_name == parent_name;
}
