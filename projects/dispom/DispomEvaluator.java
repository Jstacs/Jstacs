/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package projects.dispom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.AbstractMap.SimpleEntry;

import projects.dispom.PFMComparator.NormalizedEuclideanDistance;
import de.jstacs.classifier.ScoreBasedPerformanceMeasureDefinitions;
import de.jstacs.classifier.scoringFunctionBased.gendismix.GenDisMixClassifier;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.FileManager;
import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.motifDiscovery.MotifDiscoveryAssessment;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder.RandomSeqType;
import de.jstacs.scoringFunctions.NormalizableScoringFunction;
import de.jstacs.scoringFunctions.NormalizedScoringFunction;
import de.jstacs.scoringFunctions.directedGraphicalModels.BayesianNetworkScoringFunction;
import de.jstacs.scoringFunctions.mix.StrandScoringFunction;
import de.jstacs.scoringFunctions.mix.motifSearch.HiddenMotifsMixture;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.REnvironment;

/**
 * Test for single hidden motifs.
 * 
 * @author Jens Keilwagen
 */
public class DispomEvaluator {
	
	private static String getConsensus( AlphabetContainer con, double[][] pwm ) {
		String c = "";
		for( int m, p, l = 0; l < pwm.length; l++ ) {
			m = 0;
			for( p = 1; p < pwm[l].length; p++ ) {
				if( pwm[l][m] < pwm[l][p] ) {
					m = p;
				}
			}
			c += con.getSymbol( l, m );
		}
		return c;
	}
	
	/**
	 * @param args
	 * 0 home
	 * 1 best path
	 * 2 data path
	 * 3 ignore
	 * 4 sep
	 * 5 [optional] p-value;
	 */
	public static void main(String[] args) throws Exception {
		System.out.println( Arrays.toString( args ) );
		
		double sign = 1E-3;
		if( args.length == 6 ) {
			sign = Double.parseDouble( args[5] );
		}
		
		REnvironment r = null;
		try{
			r = new REnvironment( "localhost", "", "" );
			r.voidEval( "require( seqLogo );" );
		} catch( Exception e ) {
			System.out.println( "Could not open local Rserve-connection." );
		}
		try{	
			ArrayList<SimpleEntry<String, double[][]>> list = PFMComparator.readPFMsFromEMBL( "./transfac.dat", Integer.MAX_VALUE );
			
			File f = new File( args[0] );
			System.out.println( f.getAbsolutePath() );
			
			File[] all = f.listFiles( (FilenameFilter) new RegExFilenameFilter("xml",false,true,".*\\.xml") );
			System.out.println( "files: " + all.length );
			Arrays.sort( all, new FileComparator() );
			
			FileAssigner fa = new FileAssigner();
			double best, current;
			String n;
			char ignore = args[3].charAt( 0 );
			String[] fNames = new String[2];
			for( int idx, j = 0; j < all.length; ) {
				
				best = Double.NEGATIVE_INFINITY;
				idx = -1;
				do {
					current = infos( null, args[1], all[j], false, null, null, ignore, sign, null );
					if( current > best ) {
						best = current;
						idx = j;
					}
					j++;
				}while( j < all.length && fa.compare( all[j-1], all[j] ) == 0 );
				
				if( idx > -1 ) {
					System.out.println("---------------------------------------------------------------------------------");
					FileManager.copy( all[idx].getAbsolutePath(), args[1] + "/" + all[idx].getName() );
					n = all[idx].getName();
					n=n.substring(n.indexOf("classifier-")+11, n.indexOf( "-" + args[4] ) );
					System.out.println( n );
					n = replace( n );
					System.out.println( n );
					int x=-1;
					do {
						x = n.indexOf( "-", x+1 );
						if( x > 0 ) {
							fNames[0] = args[2] + "/" + n.substring(0, x);
							fNames[1] = args[2] + "/" + n.substring(x+1);
						}
						System.out.println( x + "\t" + fNames[0] + "\t" + fNames[1] );
					} while( x >= 0 && !( (new File( fNames[0] )).exists() && (new File( fNames[1] )).exists()) );
					if( x < 0 ) {
						System.out.println( n );
						System.out.println( "Problem: " + x );
						System.exit(1);
					} else {
						System.out.println( Arrays.toString( fNames ) );
					}
					
					//String suffix = all[idx].getName();
					//int yy = suffix.indexOf( "motifs-" );
					n = all[idx].getName();
					x = n.indexOf( "-" + args[4] );
					String e = "length-";
					int xx = n.indexOf( e, x );
					xx = n.indexOf( "-", xx+e.length() );
					n = n.substring( x, xx ).replaceAll("_", "-");
					current = infos( r, args[1], all[idx], true, list, fNames, ignore, sign, n );
					System.out.println("=================================================================================");
				}
			}
		} catch( Exception e ) {
			e.printStackTrace();
		}
		if( r != null ) {
			r.close();
		}
	}
	
	private static String replace( String n ) {
		int x = -1;
		do {
			x = n.indexOf( "[", x+1 );
			if( x >= 0 ) {
				int y = n.indexOf("]", x );
				if( y > -1 ) {
					System.out.println( n.substring(x, y) );
					n = n.substring(0,x) + n.substring(x, y).replaceAll( "-", "," ) + n.substring(y);//XXX
				}
			}
		}while( x > -1 );
		return n;
	}
	
	private static double infos( REnvironment r, String home, File f, boolean show, ArrayList<SimpleEntry<String, double[][]>> list, String[] fName, char ignore, double sign, String suffix ) throws Exception {		
		GenDisMixClassifier cl = new GenDisMixClassifier( FileManager.readFile( f ) );
		HiddenMotifsMixture md = (HiddenMotifsMixture) cl.getScoringFunction( 0 );
		if( show ) {
			System.out.println();
			System.out.println( md );
		}
		double current = cl.getLastScore();
		
		System.out.print( f.getName() + "\t" + current );
		double[][][] pwm = new double[md.getNumberOfMotifs()][][];
		for( int m = 0; m < pwm.length; m++ ) {
			NormalizableScoringFunction nsf = md.getFunction( 2*m );
			StrandScoringFunction strand;
			if( nsf instanceof NormalizedScoringFunction ) {
				strand = (StrandScoringFunction) ((NormalizedScoringFunction) nsf).getFunction();
			} else {
				strand = (StrandScoringFunction) nsf;
			}
			BayesianNetworkScoringFunction motif = (BayesianNetworkScoringFunction) strand.getFunction( 0 );
			pwm[m] = motif.getPWM();
			System.out.print( "\t" + md.getMotifLength(m) + "\t" + getConsensus( cl.getAlphabetContainer(), pwm[m] ) );
		}
		System.out.println();
		
		if( list != null ) {
			int last = fName[0].lastIndexOf( "/" );
			for( int m = 0; m < pwm.length; m++ ) {
						
				DNAAlphabet dna = new DNAAlphabet();
				
				//parameters
				System.out.println( "parameter PWM" );
				System.out.println( PFMComparator.matrixToString( pwm[m] ) );
				ComparableElement<String, Double>[] ce = PFMComparator.find( dna, pwm[m], list, new NormalizedEuclideanDistance(), 7, 2, false, 0.24 );
				for( int i = 0; i < ce.length; i++ ) {
					System.out.println( i + "\t" + ce[i].getWeight() + "\t" +ce[i].getElement() );
				}
				System.out.println("-+-+-+-+-+-+-+-+-+-+-");
	
				//prediction
				SignificantMotifOccurrencesFinder smof;
				File bg = new File( fName[1] );
				if( bg.exists() ) {
					smof = new SignificantMotifOccurrencesFinder( md, new Sample( md.getAlphabetContainer(), new SparseStringExtractor( fName[1], ignore ) ), null, sign );
				} else {
					smof = new SignificantMotifOccurrencesFinder( md, RandomSeqType.PERMUTED, true, 1000, sign );
				}
				Sample data = new Sample( md.getAlphabetContainer(), new SparseStringExtractor( fName[0], ignore ) );
				//System.out.println( data.getNumberOfElements() + "\tlength=" + data.getElementLength() );
				//System.out.println( md.getLength() );
				Sample motifs = smof.getBindingSites( data, m );
				if( motifs != null ) {
					System.out.println( "(" + motifs.getNumberOfElements() + " BS vs. " + (bg.exists()?"bg":"permuted") + ", " + sign + ")" );
					pwm[m] = PFMComparator.getPFM( motifs );
					ce = PFMComparator.find( dna, pwm[m], list, new NormalizedEuclideanDistance(), 7, 2, false,  0.24 );
					System.out.println( "binding site PFM " );
					System.out.println( PFMComparator.matrixToString( pwm[m] ) );
					for( int i = 0; i < ce.length; i++ ) {
						System.out.println( i + "\t" + ce[i].getWeight() + "\t" +ce[i].getElement() );
					}
					System.out.println("-+-+-+-+-+-+-+-+-+-+-");
					
					if( r != null ) {
						//TODO seqLogo, pos
						r.createMatrix( "pwm", pwm[m] );
						String name0 = fName[0].substring( last+1 ).replaceAll("[,._]", "-"), name1 = fName[1].substring( last+1 ).replaceAll("[,._]", "-");
						String name = name0 + /*"-" + name1 +*/ suffix + "-" + m;
						//motifs.save( new File(home+"/motifs-of-" + name0 + ".txt" ) );
						r.voidEval( "print(pwm);pwm<-apply(pwm,1,function(x){x/sum(x)});print(pwm);" );
						r.plotToPDF( "seqLogo(pwm);",8,3, home + "/seqLogo-" + name + ".pdf", true );
						pwm[m] = PFMComparator.getReverseComplement( dna, pwm[m] );
						r.createMatrix( "pwm", pwm[m] );
						r.voidEval( "print(pwm);pwm<-apply(pwm,1,function(x){x/sum(x)});print(pwm);" );
						r.plotToPDF( "seqLogo(pwm);",8,3, home + "/seqLogo-" + name + "-rc.pdf", true );
						int l = md.getLength();
						String posCmd = md.getFunction(2*m+1).toString() + "\n"
							+ "h=hist(pos,breaks=seq(0,"+l+",by=10),plot=F);\n"
							+ "plot(0,0,col=0,xlim=c(0,"+l+"),ylim=c(0,max(p,h$density)),xlab=\"Position\",ylab=\"Density\",main=\"\",axes=F,cex.lab=1.25);\n"
							+ "L="+l+";z=seq(0,L-1,by=50);at=c(z,L-1);z=z-L;z[z%%100!=0]=\"\";x=c(z,-1);axis(1,at,x,cex.axis=1.25);axis(2,cex.axis=1.25);"
							+ "plot(h,freq=F,add=T);\n"
							+ "lines(l,p,col=2);";
						IntList pos = smof.getStartPositions( 0, data, m, Integer.MAX_VALUE );
						r.createVector( "pos", pos.toArray() );
						System.out.println( home + "/position-" + name + ".pdf" );
						r.plotToPDF( posCmd,8,5, home + "/position-" + name + ".pdf", true );
					}		
/*
					// artificial data set => check predictions => auc-PR
					File pos = new File( fName[0].substring( 0,last+1 ) + "pos_" + fName[0].substring( last+1 ) );
					if( pos.exists() ) {
						System.out.println();
						Sample trueMotifs = new Sample( data.getAlphabetContainer(), new SparseStringExtractor( fName[0].substring( 0, fName[0].indexOf( "_e.txt" )+ 6 ), ignore ) );
						//System.out.println( trueMotifs );
						System.out.println( "check against motif annotation with motif length " + trueMotifs.getElementLength() );
						Sample truth = addPosition( data, pos.getAbsolutePath(), trueMotifs.getElementLength() );//TODO
						double[][] val = smof.getValuesForEachNucleotide(truth, 0, false);//TODO
						double[][] sortedVal = MotifDiscoveryAssessment.getSortedValuesForMotifAndFlanking( truth, val, smof.getOffsetForAucPR(), smof.getFactorForAucPR(), "motif " + 0 );//TODO
						double auc = ScoreBasedPerformanceMeasureDefinitions.getAUC_PR( sortedVal[0], sortedVal[1], null );
						
						System.out.println( "(nucleotide) auc-PR: " + auc );
					}
*/
				} else {
					System.out.println( "no binding sites predicted" );
				}
			}
		} 
		return current;
	}
	
	private static Sample addPosition( Sample data, String posFileName, int length ) throws Exception
	{
		BufferedReader r = new BufferedReader( new FileReader( posFileName ) );
		String line;
		Sequence[] seqs = data.getAllElements();
		for( int i = 0; i < seqs.length; i++ )
		{
			//System.out.println(i);
			line = r.readLine().trim();
			if( !line.equalsIgnoreCase( "na" ) )
			{
				seqs[i] = seqs[i].annotate( false, new MotifAnnotation( "motif 0", Integer.parseInt( line )-1, length, Strand.UNKNOWN ) );
			}
		}
		r.close();
		return new Sample( "annotated " + data.getAnnotation(), seqs );
	}
}

class FileComparator implements Comparator<File> {
	public int compare( File o1, File o2 ) {
		return o1.getName().compareTo( o2.getName() );
	}		
}

class FileAssigner implements Comparator<File> {
	static final String[] INFIX = { "best" , "enum" , "heuristic"};
	static final int getIndex( String n ) {
		int i = -1, j = 0;
		while( j < INFIX.length && (i = n.lastIndexOf( INFIX[j] )) <0 ) {
			j++;
		}
		return i;
	}
	public int compare( File o1, File o2 ) {
		String n1 = o1.getName();
		String n2 = o2.getName();
		//System.out.println( n1 + "\t" + n2 );
		try {
			n1 = n1.substring( 0, getIndex(n1) );
			n2 = n2.substring( 0, getIndex(n2) );
		} catch( RuntimeException e ) {
			System.err.println( n1 );
			System.err.println( n2 );
			throw e;
		}
		//System.out.println( n1 + "\t" + n2 );
		return n1.compareTo( n2 );
	}		
}