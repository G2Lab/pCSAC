#ifndef XYPDBATOM_H
#define XYPDBATOM_H

#include <vector>
#include <algorithm>

//#include "XYUtility.h"
#include "XYPoint3D.h"
/*
	http://www.wwpdb.org/documentation/format23/sect9.html#ATOM
	ATOM
	Overview

	The ATOM records present the atomic coordinates for standard residues (see http://deposit.pdb.org/public-component-erf.cif). They also present the occupancy and temperature factor for each atom. Heterogen coordinates use the HETATM record type. The element symbol is always present on each ATOM record; segment identifier and charge are optional.
	Record Format

	COLUMNS      DATA TYPE        FIELD      DEFINITION
	------------------------------------------------------
	 1 -  6      Record name      "ATOM    "
	 7 - 11      Integer          serial     Atom serial number.
	13 - 16      Atom             name       Atom name.
	17           Character        altLoc     Alternate location indicator.
	18 - 20      Residue name     resName    Residue name.
	22           Character        chainID    Chain identifier.
	23 - 26      Integer          resSeq     Residue sequence number.
	27           AChar            iCode      Code for insertion of residues.
	31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
											 Angstroms
	39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
											 Angstroms
	47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
											 Angstroms
	55 - 60      Real(6.2)        occupancy  Occupancy.
	61 - 66      Real(6.2)        tempFactor Temperature factor.
	77 - 78      LString(2)       element    Element symbol, right-justified.
	79 - 80      LString(2)       charge     Charge on the atom.

	Details

	# ATOM records for proteins are listed from amino to carboxyl terminus.
	# Nucleic acid residues are listed from the 5' to the 3' terminus.
	# No ordering is specified for polysaccharides.
	# The list of ATOM records in a chain is terminated by a TER record.
	# If more than one model is present in the entry, each model is delimited by MODEL and ENDMDL records.
	# If an atom is provided in more than one position, then a non-blank alternate location indicator must be used as the alternate location indicator for each of the positions. Within a residue all atoms that are associated with each other in a given conformation are assigned the same alternate position indicator.
	# For atoms that are in alternate sites indicated by the alternate site indicator, sorting of atoms in the ATOM/ HETATM list uses the following general rules:

	In the simple case that involves a few atoms or a few residues with alternate sites, the coordinates occur one after the other in the entry.

	In the case of a large heterogen groups which are disordered, the atoms for each conformer are listed together.

	# The insertion code is commonly used in sequence numbering
	# If the depositor provides the data, then the isotropic B value is given for the temperature factor.
	# If there are neither isotropic B values from the depositor, nor anisotropic temperature factors in ANISOU, then the default value of 0.0 is used for the temperature factor.
	# Columns 77 - 78 contain the atom's element symbol (as given in the periodic table), right-justified.
	# Columns 79 - 80 indicate any charge on the atom, e.g., 2+, 1-. In most cases these are blank.

	Verification/Validation/Value Authority Control

	PDB checks ATOM/HETATM records for PDB format, sequence information, and packing. The PDB reserves the right to return deposited coordinates to the author for transformation into PDB format.
	Relationships to Other Record Types

	The ATOM records are compared to the corresponding sequence database. Residue discrepancies appear in the SEQADV record. Missing atoms are annotated in the remarks. HETATM records are formatted in the same way as ATOM records. The sequence implied by ATOM records must be identical to that given in SEQRES, with the exception that residues that have no coordinates, e.g., due to disorder, must appear in SEQRES.
	Example

			 1         2         3         4         5         6         7         8
	12345678901234567890123456789012345678901234567890123456789012345678901234567890
	ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92           N
	ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85           C
	ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34           C
	ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65           O
	ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88           C
	ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41           C
	ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64           C
	ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11           C
	ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58           C
	ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25           C

	Known Problems

	No distinction is made between ribo- and deoxyribonucleotides in the SEQRES records. These residues are identified with the same residue name (i.e., A, C, G, T, U). 

	Another link should be useful.
	# http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html#pdb-atom-name-anomalies
	Notes on atom naming
	The atom name is a four-character field that may be subdivided into atomic symbol (A2), "remoteness indicator" (A1), and "branch designator" (A1) subfields, as detailed below. The atom name may be subdivided into the following subfields:

	Name 	Format 	Description
	atomic symbol 	A2 	Right-justified atomic symbol, e.g. " C".
	remoteness indicator 	A1 	Greek letter distance abbreviation. In order of increasing distance, these are "A" for alpha, "B" for beta, "G" for gamma, "D" for delta, "E" for epsilon, "Z" for zeta, and "H" for eta.
	branch designator 	A1 	digit designating the branch direction; left blank if the sidechain is unbranched.

	For the name " CG1;", " C" denotes the species (carbon), "G" identifies it as a gamma atom, and "1" denotes the branch of a beta-branched amino acid. Note that the atomic symbol is right-justified (i.e. " C", "ZN"), so " CA " is an alpha carbon, and "CA  " is a calcium atom.

	For amino acids, the traditional order of atom records within a residue is N, CA, C, O, followed by the sidechain atoms (CB, CG1, CG2 . . . ) in order first of increasing remoteness, and then branch. The extra oxygen at the carboxyl terminal is called " OXT", and appears after all sidechain atoms. 
	*/
using namespace std;
class CXYPDBAtom
{
public:
	CXYPDBAtom(void);
	~CXYPDBAtom(void);
	void Initialize(void);
	void ReadLine(const char *acBuffLine);
	
	char *GetRecName();
	void SetRecName(char *acRecName);
	int	GetSerial();
	void SetSerial(int iSerial);

	char *GetName();
	void SetName(char *acAtom);
	char *GetAtomicSym();
	void SetAtomicSym(char *acAtomicSym);
	char *GetRemoteInd();
	void SetRemoteInd(char *acRemoteInd);
	char *GetBranchDes();
	void SetBranchDes(char *acBranchDes);

	char *GetAltLoc();
	void SetAltLoc(char *acAltLoc);
	char *GetResName();
	void SetResName(char *acResName);
	char *GetChainID();
	void SetChainID(char *acChainID);
	int GetResSeq();
	void SetResSeq(int iResSeq);
	char *GetICode();
	void SetICode(char *acICode);

	float GetX();
	void SetX(float fX);
	float GetY();
	void SetY(float fY);
	float GetZ();
	void SetZ(float fZ);

	float GetOccupancy();
	void SetOccupancy(float fOccupancy);
	float GetTempFactor();
	void SetTempFactor(float fTempFactor);

	char *GetSegID();
	void SetSegID(char *acSegID);

	char *GetElement();
	void SetElement(char *acElement);
	char *GetCharge();
	void SetCharge(char *acCharge);

	CXYPoint3D<float> GetPoint() const;
	void SetPoint(const CXYPoint3Df& rP);




private:
	char	m_acRecName[7];	// 1 -  6      Record name      "ATOM    "
	int		m_iSerial;		// 7 - 11      Integer          serial     Atom serial number.

	char	m_acAtom[5];		// 13 - 16      Atom             name       Atom name.
	// more detail see http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html#pdb-atom-name-anomalies
	// subdetail of Atom Name
	char	m_acAtomicSym[3];	// 13-14  
	char	m_acRemoteInd[2]; // 15
	char	m_acBranchDes[2]; // 16

	char	m_acAltLoc[2];	// 17           Character        altLoc     Alternate location indicator.
	char	m_acResName[4];	// 18 - 20      Residue name     resName    Residue name.
	char	m_acChainID[2];	// 22           Character        chainID    Chain identifier.
	int		m_iResSeq;		// 23 - 26      Integer          resSeq     Residue sequence number.
	char	m_acICode[2];		// 27           AChar            iCode      Code for insertion of residues.

	float	m_fX;			// 31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angstroms
	float	m_fY;			// 39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angstroms
	float	m_fZ;			// 47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angstroms

	float	m_fOccupancy;	// 55 - 60      Real(6.2)        occupancy  Occupancy.
	float	m_fTempFactor;	// 61 - 66      Real(6.2)        tempFactor Temperature factor.
	
	char	m_acSegID[5];	// 73 - 76		LString(4)		segment identifier, left-justified. [format version 2.0 and later.]
	
	char	m_acElement[3];	// 77 - 78      LString(2)       element    Element symbol, right-justified.
	char	m_acCharge[3];	// 79 - 80      LString(2)       charge     Charge on the atom.
	CXYPoint3Df m_point;
};


#endif



