package projects.dimont.hts;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import javax.naming.OperationNotSupportedException;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.utils.Pair;
import de.jstacs.utils.ToolBox;
import projects.dimont.AbstractSingleMotifChIPper;

public enum HTS_InitialPwms {
	
	DIMONT,
	INI_MOTIF,
	CYCLE_CORRECTION;
	
	static void compute(AbstractSingleMotifChIPper model, int motifLength, Sequence kmer, DataSet data, double[] weights,
			Pair<DataSet,double[][]> alternativeData, double pwmCorrectionFactor, double h, HTS_InitialPwms type) 
			throws IllegalArgumentException, EmptyDataSetException, WrongAlphabetException, Exception {
		double[] parameters;
		
		
		switch (type) {
		case DIMONT:
			//model.initializeMotif( 0, new DataSet("", kmer ), new double[]{h} );
			DimontTool.setMotifParameters(model, kmer, h);
			
			break;
		case INI_MOTIF:
			double[][] pwm = getPwmFromKmer(kmer.toString(), data, weights, motifLength);
			parameters = pwmToInputParameters(pwm);
			model.reset();
			model.initializeHiddenUniformly();
			model.setParametersForFunction(0, parameters, 0);
			break;
		case CYCLE_CORRECTION:
			DataSet seqs = alternativeData.getFirstElement();
			double[] weights_fg = alternativeData.getSecondElement()[0];
						
			Collection<Sequence> fg = new ArrayList<Sequence>();
			Collection<Sequence> bg = new ArrayList<Sequence>();
			for (int i=0;i<seqs.getNumberOfElements() ;i++) {
				if (weights_fg[i]==1.0) {
					fg.add(seqs.getElementAt(i));
				}
				else {
					bg.add(seqs.getElementAt(i));
				}
			}
			
			double[] w = new double[fg.size()];
			Arrays.fill(w, 1.0);
			double[][] pwm_fg = getPwmFromKmer(kmer.toString(), new DataSet("fg", fg), w, motifLength);
			w = new double[bg.size()];
			Arrays.fill(w, 1.0);
			double[][] pwm_bg = getPwmFromKmer(kmer.toString(), new DataSet("bg", bg), w, motifLength);
			
			// pwm_fg = pwm_fg - correctionFactor*pwm_bg
			for (int col=0;col<pwm_fg.length;col++) {
				for (int row=0;row<pwm_fg[0].length;row++) {
					pwm_fg[col][row] = pwm_fg[col][row] - pwmCorrectionFactor*pwm_bg[col][row];
					if (pwm_fg[col][row]<=0) pwm_fg[col][row] = 0.01;
				}
				double sum = ToolBox.sum(pwm_fg[col]);
				for (int row=0;row<pwm_fg[0].length;row++) {
					pwm_fg[col][row] /= sum;

				}
			}
			
//			System.out.println("\n");
//			for (int row=0;row<pwm_fg[0].length;row++) {
//				String roww =""; 
//				for (int col=0;col<pwm_fg.length;col++) {
//					roww = roww + "\t" + pwm_fg[col][row];
//				}
//				System.out.println(roww);
//			}		
			
			parameters = pwmToInputParameters(pwm_fg);
			model.reset();
			model.initializeHiddenUniformly();
			model.setParametersForFunction(0, parameters, 0);
			break;
		}

	}
	
	private static double[] pwmToInputParameters(double[][] pwm) {
		double[] parameters = new double[pwm.length*pwm[0].length];
		for (int col=0;col<pwm.length;col++) {
			for (int row=0;row<pwm[0].length;row++) {
				parameters[col*pwm[0].length+row] = Math.log(pwm[col][row]);
			}
		}
		return parameters;
	}
	
	private static double[][] getPwmFromKmer(String kmer, DataSet data, double[] weights, int motifLength) throws WrongAlphabetException, OperationNotSupportedException{
		AlphabetContainer con = data.getAlphabetContainer();
		if( !con.isSimple() || !con.isDiscrete() ) {
			throw new WrongAlphabetException();
		}
		
		DiscreteAlphabet alphabet = (DiscreteAlphabet) data.getAlphabetContainer().getAlphabetAt(0);
		Hashtable<String, String> relatedKmers = new Hashtable<String, String>(2*kmer.length()*((int) alphabet.length()-1)+2,1);
		Hashtable<String, Double> weightedKmerCounts = new Hashtable<String, Double>(kmer.length()*((int) alphabet.length()-1)+1,1);
		
		// define kmers with Hamming distance of 1
		String[][] sequences = new String[kmer.length()][(int) alphabet.length()];
		for (int col=0; col<kmer.length(); col++) {
			for (int row=0; row<alphabet.length(); row++) {
				String related = kmer.substring(0, col)+alphabet.getSymbolAt(row)+kmer.substring(col+1, kmer.length());
				String revComplSequence =  Sequence.create(con, related).reverseComplement().toString();
				sequences[col][row] = related;
				relatedKmers.put(related, related);
				relatedKmers.put(revComplSequence, related);
				weightedKmerCounts.put(related, 0.0);
			}
		}
		
		int k = kmer.length();
		boolean kmerFound = false;
		String newKmer = "", lastKmer = "";
		int lastStart, startPosition;
		Sequence seq;
		String s;
		// count defined kmers (only sequences that contain exactly one defined kmer)
		// run over all sequences
		for (int n=0; n<weights.length; n++) {
			seq = data.getElementAt( n );
			s = seq.toString();
			lastStart = seq.getLength()-k;
			
			//run over all k-mers
			kmerFound = false;
			lastKmer = "";
			for( startPosition = 0; startPosition <= lastStart; startPosition++ ) {
				newKmer = s.substring( startPosition, startPosition+k );
				if( relatedKmers.containsKey( newKmer ) ) {
					if (lastKmer.equals("") || lastKmer.equals(relatedKmers.get(newKmer))) {
						kmerFound = true;
						lastKmer = relatedKmers.get(newKmer);
						startPosition +=2;
					}
					else {
						kmerFound = false;
						break;
					}
				}
			}
			if (kmerFound) {
				weightedKmerCounts.put( lastKmer, weightedKmerCounts.get(lastKmer)+weights[n] );
			}
		}
		
		// create pwm
		double[][] pwm = new double[motifLength][(int) alphabet.length()];
		int numBorderingPositions = motifLength - k;
		int motifStart = (int) Math.floor(numBorderingPositions / 2.0);
		// motif
		for (int col=motifStart; col<motifStart+kmer.length(); col++) {
			for (int row=0; row<alphabet.length(); row++) {
				pwm[col][row] = weightedKmerCounts.get(sequences[col-motifStart][row]) +1.0; // pseudocount 1
			}
			double sum = ToolBox.sum(pwm[col]);
			for (int row=0; row<alphabet.length(); row++) {
				pwm[col][row] /= sum;
			}
		}
		// left margin
		for(int col=0; col<motifStart; col++) {
			for (int row=0; row<alphabet.length(); row++) {
				pwm[col][row] = 0.25;
			}
		}
		// right margin
		for(int col=motifStart+kmer.length(); col<motifLength; col++) {
			for (int row=0; row<alphabet.length(); row++) {
				pwm[col][row] = 0.25;
			}
		}
		
		return pwm;
	}
	
	

}
