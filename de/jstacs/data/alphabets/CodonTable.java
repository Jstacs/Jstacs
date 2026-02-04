package de.jstacs.data.alphabets;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;

/**
 * This class implements a codon table, which allows to map from codons to amino acids.
 * 
 * @author Jens Keilwagen
 */
public class CodonTable {

	String[] aminoAcid;
	int[][][] mapping;
	
	/**
	 * Constructor for the standard genetic code.
	 * 
	 * @throws IOException if something went wrong while processing the file
	 * @throws WrongAlphabetException if the used alphabet is not compatible with DNA
	 * 
	 * @see DNAAlphabet#SINGLETON
	 */
	public CodonTable() throws IOException, WrongAlphabetException {
		this("projects/gemoma/test_data/genetic_code.txt");
	}

	/**
	 * General constructor allowing to used arbitrary codon mappings given in a file.
	 * 
	 * @param fName the file name
	 * 
	 * @throws IOException if something went wrong while processing the file
	 * @throws WrongAlphabetException if the used alphabet is not compatible with DNA
	 * 
	 * @see DNAAlphabet#SINGLETON
	 */
	public CodonTable( String fName ) throws IOException, WrongAlphabetException {
		//prepare
		ArrayList<String> aa = new ArrayList<String>();
		mapping = new int[4][4][4];
		for( int i = 0; i < mapping.length; i++ ) {
			for( int j = 0; j < mapping[i].length; j++ ) {
				Arrays.fill(mapping[i][j], -1);
			}
		}
		
		//read from file
		DNAAlphabet dna = DNAAlphabet.SINGLETON;
		BufferedReader r = new BufferedReader( new FileReader(fName) );
		String line;
		while( (line=r.readLine()) != null ) {
			String[] s = line.split("\t");
			if( aa.contains(s[0]) ) {
				throw new IllegalArgumentException( "Amino acid "+ s[0] +" was used multiple times." );				
			}
			aa.add(s[0]);
			String[] triplets = s[1].split(",");
			for( String t: triplets ) {
				t = t.trim();
				int i = dna.getCode(t.substring(0, 1));
				int j = dna.getCode(t.substring(1, 2));
				int k = dna.getCode(t.substring(2, 3));
				if( mapping[i][j][k]!=-1 ) {
					//already used
					throw new IllegalArgumentException( "Codon "+ t +" was used multiple times." );
				} else {
					mapping[i][j][k] = aa.size()-1;
				}
			}
		}
		r.close();
		aminoAcid = aa.toArray(new String[0]);
		
		//check
		for( int i = 0; i < mapping.length; i++ ) {
			for( int j = 0; j < mapping[i].length; j++ ) {
				for( int k = 0; k < mapping[j].length; k++ ) {
					if( mapping[i][j][k] == -1 ) {
						throw new IllegalArgumentException("mapping not complete");
					}
				}
			}
		}
	}
	
	/**
	 * Return the number of different amino acids (including stop).
	 * 
	 * @return the number of different amino acids
	 */
	public int getNumberOfAminoAcids() {
		return aminoAcid.length;
	}

	/**
	 * Returns the String representations of all amino acids that are used.
	 * This method might be used to create a {@link DiscreteAlphabet}. 
	 * 
	 * @return the String representations of all amino acids that are used.
	 * 
	 * @see DiscreteAlphabet#DiscreteAlphabet(boolean, String...)
	 */
	public String[] getAminoAcids() {
		return aminoAcid.clone();
	}
	
	/**
	 * Return the String representation of the amino acid with a given index.
	 * 
	 * @param idx the index of the amino acid
	 * 
	 * @return the String representation of the amino acid with a given index
	 */
	public String getRepresentation( int idx ) {
		return aminoAcid[idx];
	}
	
	/**
	 * Returns the index of the amino acid encoded in the Sequence <code>seq</code> starting at position <code>start</code>.

	 * @param seq the sequence
	 * @param start the start position
	 * 
	 * @return the index of the encoded amino acid
	 */
	public int getIndex( Sequence seq, int start ) {
		//if( !seq.getAlphabetContainer().checkConsistency( DNAAlphabetContainer.SINGLETON ) ) throw new WrongAlphabetException();
		return mapping[seq.discreteVal(start)][seq.discreteVal(start+1)][seq.discreteVal(start+2)];
	}

	/**
	 * Returns the String representation of the amino acid encoded in the Sequence <code>seq</code> starting at position <code>start</code>.

	 * @param seq the sequence
	 * @param start the start position
	 * 
	 * @return the String of the encoded amino acid
	 * 
	 * @see #getIndex(Sequence, int)
	 * @see #getRepresentation(int)
	 */
	public String getRepresentation( Sequence seq, int start ) {
		return getRepresentation( getIndex(seq, start) );
	}
}