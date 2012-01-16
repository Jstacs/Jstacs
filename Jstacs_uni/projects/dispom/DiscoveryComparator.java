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
import java.io.FileFilter;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

import javax.naming.OperationNotSupportedException;

import projects.dispom.PFMComparator.PFMDistance;
import de.jstacs.NonParsableException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.classifier.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifier.differentiableSequenceScoreBased.ScoreClassifier;
import de.jstacs.classifier.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.classifier.differentiableSequenceScoreBased.msp.MSPClassifier;
import de.jstacs.classifier.performanceMeasures.PRCurve;
import de.jstacs.classifier.performanceMeasures.ROCCurve;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.annotation.MotifAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.StrandedLocatedSequenceAnnotationWithLength.Strand;
import de.jstacs.io.AbstractStringExtractor;
import de.jstacs.io.FileManager;
import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.motifDiscovery.MotifDiscoveryAssessment;
import de.jstacs.motifDiscovery.SignificantMotifOccurrencesFinder;
import de.jstacs.results.ResultSet;
import de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels.structureLearning.measures.Measure;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.DurationDiffSM;
import de.jstacs.sequenceScores.statisticalModels.differentiable.mix.motif.ExtendedZOOPSDiffSM;
import de.jstacs.utils.ComparableElement;
import de.jstacs.utils.IntList;
import de.jstacs.utils.REnvironment;

/** 
 * A comparator for the results of different de-novo motif discoverer.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class DiscoveryComparator {

	private static final String HOME = "../Dispom/";//"../Thesis/results/Dispom/";//TODO
		// subfolder: data and results
	private static final String TOOLS = "tools/";
	private static final String DATATYPE = "artificial/";
	private static final String motifName = "implanted motif";
	private static final String host = "localhost";
	private static final String user = "";
	private static final String pass = "";
	
	private static enum Discoverer {
		A__GLAM("black"),
		//AMADEUS("lawngreen"),
		DEME("blue"),
		DME("cyan"),
		//FIRE,
		GIBBS_SAMPLER("orange"),
		//CENTROID_GIBBS_SAMPLER(""),
		IMPROBIZER("skyblue1"),
		MEME("red"),
		MOAN("yellow"),
		WEEDER("magenta2"),
		Dispom("green3");
		
		private static int no = 0;
		private int myNo, used;
		private String color;
		
		Discoverer( String color ) {
			setMyNo();
			used = 3;
			this.color = color;
		}
		
		private void setMyNo() {
			myNo = 5 + Discoverer.no++;
		}
		
		public int getNumber() {
			return myNo;
		}
		
		public void used() {
			used++;
		}
		
		public static void resetUsed(){
			Discoverer[] vals = Discoverer.values();
			for(Discoverer d : vals){
				d.used=3;
			}
		}
		
		public String getColor() {
			return color;
		}
	}
	
	private static class DiscoveryResult {
		private static PRCurve prCurve = new PRCurve();
		private static ROCCurve rocCurve = new ROCCurve();
		
		ResultSet rs;
		DoubleTableResult pr;
		DoubleTableResult roc;
		String description;
		Discoverer d;
		double[][] pwm;
		ExtendedZOOPSDiffSM md;
		DataSet truth;
		double aucPR;
		int pch;
		
		private DiscoveryResult( Discoverer d, DataSet truth, DataSet data, boolean remove, String... strings ) throws Exception {
			this.d = d;
			pch = d.used;
			d.used();
			this.truth = truth;
			DataSet pred = null;
			SignificantMotifOccurrencesFinder smof = null;
			//load complete predictions
			switch( d ) {
				case A__GLAM:
					pred = ParserToolBox.annotateWithAGlamResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "aglam/" + strings[0], "motif 1",true );	
					break;
				/*
				case AMADEUS:
					pred = ParserToolBox.annotateWithAmadeusResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "amadeus/" + strings[0], 1, Integer.parseInt( strings[1] ) );
					break;
				*/
				case DEME:
				    pred = ParserToolBox.annotateWithDemeResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "deme/" + strings[0] );
				    break;
				case DME:
					pred = ParserToolBox.annotateWithDMEResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "dme/" + strings[0], 0 );
					break;
				case GIBBS_SAMPLER:
					pred = ParserToolBox.annotateWithGibbsSamplerResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "gibbs/" + strings[0] );
					break;
				/*	
				case CENTROID_GIBBS_SAMPLER:
					//XXX
					pred = ParserToolBox.annotateWithGibbsSamplerResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "centroid-gibbs/" + strings[0] );
					break;
				*/
				case IMPROBIZER:    
				    pred = ParserToolBox.annotateWithImprobizerResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "improbizer/" + strings[0], ParserToolBox.getImprobizerMotifLength( HOME + "results/" + DATATYPE + TOOLS + "improbizer/" + strings[1] ) );
				    break;
				case MEME:
					pred = ParserToolBox.annotateWithMemeResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "meme/" + strings[0] );
					break;
				case MOAN:
					pred = ParserToolBox.annotateWithMoAnResults( data, false, (new File(HOME + "results/" + DATATYPE + TOOLS + "moan/")).listFiles( (FileFilter) new RegExFilenameFilter( "", false, true, strings[0] ) ) );
					break;
				case WEEDER:
				    pred = ParserToolBox.annotateWithWeederResults( data, false, HOME + "results/" + DATATYPE + TOOLS + "weeder/" + strings[0] /*, Integer.parseInt( strings[1] )*/ );//TODO
				    break;
				case Dispom:
					ScoreClassifier cl;
					//System.out.println(HOME + "data/" + DATATYPE +prefix+".sites_"+org+".txt_e.txt_"+inf+".txt -> "+f.file.getAbsolutePath());
						
					//wr.println( f.getName() );
					File f = new File( strings[0] );
					try{
						cl = new MSPClassifier( FileManager.readFile( f ) );
					}catch(NonParsableException e){
						cl = new GenDisMixClassifier( FileManager.readFile( f ) );
					}
					md = (ExtendedZOOPSDiffSM) cl.getScoringFunction( 0 );

					smof = new SignificantMotifOccurrencesFinder( md, new DNADataSet(strings[1],AbstractStringExtractor.USUALLY), null, Double.parseDouble( strings[2] ) );
					//smof = new SignificantMotifOccurrencesFinder(md,RandomSeqType.PERMUTED,1000,sign);
					
					pred = smof.annotateMotif( data, 0 );				
					break;
				default:
					System.out.println( "not found" );
					System.exit( 1 );
			}
			if( d == Discoverer.IMPROBIZER && d == Discoverer.WEEDER /*&& d == Discoverer.AMADEUS*/ ) {
				description = strings[0] + " and " + strings[1];
			} else {
				description = strings[0];
			}
			rs = MotifDiscoveryAssessment.assess( truth, pred, 5 ).getValue()[0];
			
			if(rs.getNumberOfResults() > 0){
				//curves
				double[][] scores;
				if( d == Discoverer.Dispom ) {
					scores = MotifDiscoveryAssessment.getSortedValuesForMotifAndFlanking( truth,
							smof.getValuesForEachNucleotide( data, 0, false ),
							smof.getOffsetForAucPR(), smof.getFactorForAucPR(), motifName );
				} else {
					scores = MotifDiscoveryAssessment.getSortedScoresForMotifAndFlanking( truth, pred, motifName );
				}
				
				LinkedList<double[]> listOfDoubles = new LinkedList<double[]>();
				
				ResultSet curve = prCurve.compute( scores[0], scores[1] ); 
				aucPR = (Double) curve.getResultForName( "AUC-PR" ).getValue();
				System.out.println( getName() + "\t" + aucPR );
				pr = (DoubleTableResult) curve.getResultForName( "PR curve" ).getValue();
				if( remove ) {
					double[][] temp = pr.getValue();
					for( int i = 0; i < temp.length; i++ ) {
						if(temp[i][0] <= (Double) rs.getResultAt( 0 ).getValue()){
							listOfDoubles.add( temp[i] );
						}
					}
					pr = new DoubleTableResult( pr.getName(), pr.getComment(), listOfDoubles );
				}
			
				curve = rocCurve.compute( scores[0], scores[1] );
				roc = (DoubleTableResult) curve.getResultForName( "ROC curve" ).getValue();
				if( remove ) {
					listOfDoubles.clear();
					double[][] temp = pr.getValue();
					for( int i = 0; i < temp.length; i++ ) {
						if(temp[i][0] <= (Double) rs.getResultAt( 0 ).getValue()){
							listOfDoubles.add( temp[i] );
						}
					}
					roc = new DoubleTableResult( roc.getName(), roc.getComment(), listOfDoubles );
				}
				
				/*
				//filter strong binding sites
				switch( d ) {
					case A__GLAM:
					case AMADEUS:
					case DEME:
					case DME:
					case GIBBS_SAMPLER:
					case CENTROID_GIBBS_SAMPLER:
					case IMPROBIZER:    
					case MEME:
						break;
					case WEEDER:
						ArrayList<Sequence> list = new ArrayList<Sequence>();
						Sequence seq;
						ArrayList<SequenceAnnotation> annot = new ArrayList<SequenceAnnotation>();
						double s;
						for( int i = 0; i < pred.getNumberOfElements(); i++ ) {
							seq = pred.getElementAt( i );
							SequenceAnnotation[] seqAn = seq.getAnnotation();
							annot.clear();
							for( int j = 0; j < seqAn.length; j++ ) {
								s = (Double) (seqAn[j].getResultAt( 0 )).getResult();
								if( s >= 85 ) {
									annot.add( seqAn[j] );
								}
							}
							if( annot.size() > 0 ) {
								list.add( seq.annotate( false, annot.toArray( new SequenceAnnotation[0] ) ) );
							}
						}
						if( list.size() == 0 ) {
							pred = null;
						} else {
							pred = new Sample( "filtered data set", list.toArray( new Sequence[0] ) );
						}
					    break;
					case DiPoMM:
						break;
					default:
						System.out.println( "not found" );
				}
				/**/
			}else{
				System.out.println("NO RESULT FOR "+d.toString() + " (" + strings[0] +")" );
				pred = null;
				aucPR = 0;
			}
			
			//compute PWM
			if( pred == null ) {
				pwm = new double[0][0];
			} else {
				pwm = DiscoveryComparator.getPWM( pred, null, false );
			}
		}
		
		public double[][] getPWM(){
			return pwm;
		}
		
		public double getPPVFor(double Sn){
			if( pr == null || rs.getNumberOfResults() == 0 ) {
				return Double.NaN;
			} else {				
				int i=pr.getNumberOfLines()-1;
				while( i >= 0 && pr.getLine( i )[0] < Sn ){
					i--;
				}
				if(i == -1){
					return Double.NaN;
				}else{
					return pr.getLine( i )[1];
				}
			}
		}
		
		public boolean plot() {
			if( pr == null ) {
				return false;
			} else {
				int i=pr.getNumberOfLines()-1;
				if( pr.getLine(i)[0] > 0.1 ) {
					return true;
				} else {
					while( i >= 0 && pr.getLine( i )[1] < 0.1 ){
						i--;
					}
					return i != -1;
				}
			}
		}
		
		public String getColor() {
			return d.getColor();
		}
		
		public int getPCH() {
			return pch;
		}		
		
		public String getName() {
			String s = d.name();
			if( d == Discoverer.Dispom && pch != 3 ) {
				s += "-" + (pch-3);
			}
			return s.replaceAll( "__", "-" ).replace( '_', ' ' );
		}
		
		public double getAUC_PR() {
			return aucPR;
		}
		
		public double getValue( int index ) {
			return (Double) rs.getResultAt( index ).getValue();
		}
		
		public DoubleTableResult getPrCurve(){
			return pr;
		}
		
		//seqLogo, ...
		public void createPlots( REnvironment r, boolean rc, String suffix ) throws Exception {
			plotSeqLogo( r, HOME+"results/" + DATATYPE + "seqLogos/seqLogo-"+getName().replaceAll( "[\\._\\s]", "-" )+"-"+suffix+".pdf", rc? pwm : PFMComparator.getReverseComplement(DNAAlphabet.SINGLETON, pwm) );
			if( d == Discoverer.Dispom ) {
				plotPositionalDistribution( r, HOME+"results/" + DATATYPE + "positionalDistribution/pos-"+getName().replaceAll( "[\\._\\s]", "-" )+"-"+suffix+".pdf", md, truth );
			}
		}
		
		public DoubleTableResult getROCCurve(){
			return roc;
		}
		
		public String toString() {
			StringBuffer sb = new StringBuffer(300);
			sb.append( d.name() + "\t" );
			for( int i = 0; i < rs.getNumberOfResults(); i++ ) {
				sb.append( rs.getResultAt(i).getValue() + "\t" );
			}
			sb.append( description );
			return sb.toString();
		}
		
		public String getHeading(){
			StringBuffer sb = new StringBuffer(300);
			sb.append( "discoverer\t" );
			for( int i = 0; i < rs.getNumberOfResults(); i++ ) {
				sb.append( rs.getResultAt(i).getName() + "\t" );
			}
			sb.append( "description" );
			return sb.toString();
		}
	}
	
	private static DataSet addPosition( DataSet data, String posFileName, String revcomFileName, int length ) throws Exception
	{
		BufferedReader r = new BufferedReader( new FileReader( posFileName ) );
		BufferedReader rr = new BufferedReader( new FileReader( revcomFileName ) );
		String line,liner;
		Sequence[] seqs = data.getAllElements();
		for( int i = 0; i < seqs.length; i++ )
		{
			//System.out.println(i);
			line = r.readLine().trim();
			liner = rr.readLine().trim();
			if( !line.equalsIgnoreCase( "na" ) )
			{
				Strand strand = liner.equalsIgnoreCase( "true" ) ? Strand.REVERSE : Strand.FORWARD;
				seqs[i] = seqs[i].annotate( false, new MotifAnnotation( motifName, Integer.parseInt( line )-1, length, strand ) );
			}
		}
		r.close();
		rr.close();
		return new DataSet( "annotated " + data.getAnnotation(), seqs );
	}
	
	private static ComparableElement<Boolean, Double> getMeanDivergence( double[] bg, double[][] pwm, double[][] pwmRefFw, double[][] pwmRefRev, PFMDistance dist ){
		double forward = getMeanDivergence( bg, pwm, pwmRefFw, dist );
		double reverse = getMeanDivergence( bg, pwm, pwmRefRev, dist );
		if( forward < reverse ) {
			return new ComparableElement<Boolean, Double>( true, forward );
		} else {
			return new ComparableElement<Boolean, Double>( false, reverse );
		}
	}
	
	private static double getMeanDivergence( double[] bg, double[][] pwm, double[][] pwmRef, PFMDistance dist ){
		if(pwm.length < pwmRef.length){
			pwm = fillPwm(pwm,pwmRef.length,bg);
		}else if(pwm.length > pwmRef.length){
			pwmRef = fillPwm(pwmRef,pwm.length,bg);
		}
		return dist.compare( pwm, pwmRef, pwmRef.length/2 );
	}
	
	private static double[][] fillPwm(double[][] original, int newLength, double[] bg){
		double[][] newPwm = new double[newLength][];
		int i=0;
		for(;i<(newLength-original.length)/2;i++){
			newPwm[i] = bg.clone();
		}
		for(int k=0;k<original.length;k++,i++){
			newPwm[i] = original[k].clone();
		}
		for(;i<newLength;i++){
			newPwm[i] = bg.clone();
		}
		return newPwm;
	}
	
	private static void plotPositionalDistribution( REnvironment r, String file, ExtendedZOOPSDiffSM mix, DataSet truth) throws Exception{
		r.createVector( "pos", getPositions( truth, motifName ) );
		DurationDiffSM dsf = (DurationDiffSM) mix.getFunction( 1 );
		String plotcmd = dsf.toString()+"\n" +
				"par(mar=c(5,6,0,0),cex.lab=2,cex.axis=2);\n" +
				"hist(pos,xlab=\"Position\", ylab=\"\",main=\"\", xlim=c(min(l),max(l)),freq=F,breaks=seq(0,500,by=20),col=gray(0.8),axes=F);\n" +
				"a=axis(1,col.axis=0);b=a-501;axis(1,at=a,label=b);axis(2,las=2);\n" +
				"lines(l,p,lwd=2);";
		//System.out.println(plotcmd);
		r.plotToPDF( plotcmd, 8, 3, file, true );
	}
	
	private static void plotSeqLogo( REnvironment r,String file, double[][] pwm) throws Exception{
		if(pwm.length == 0){
			return;
		}
		r.createMatrix( "pwm", pwm );
		r.voidEval( "print(pwm);pwm<-apply(pwm,1,function(x){x/sum(x)});print(pwm);" );
		r.plotToPDF( "library(seqLogo);seqLogo(pwm);",8,3, file, true );
	}
	
	private static double[][] getPWM( DataSet data, String motifName, boolean rc ) throws Exception{
		try{
			DataSet ex = exciseSites( data, motifName, rc );

			return PFMComparator.getPFM( ex );
			
		}catch(EmptyDataSetException ex){
			return new double[0][0];
		}
	}
	
	private static int[] getPositions( DataSet data, String motifName ){
		IntList list = new IntList();
		
		for(Sequence seq : data){
			SequenceAnnotation[] anns = seq.getAnnotation();
			for(int i=0;anns != null && i<anns.length;i++){
				if(anns[i] instanceof MotifAnnotation && (motifName == null || motifName.equals( anns[i].getIdentifier() ) ) ){
					list.add( ((MotifAnnotation)anns[i]).getPosition() );
				}
			}
		}
		
		return list.toArray();
	}
	
	private static DataSet exciseSites( DataSet data, String motifName, boolean rc ) throws EmptyDataSetException, WrongAlphabetException, OperationNotSupportedException{
		LinkedList<Sequence> seqs = new LinkedList<Sequence>();
		
		for(Sequence seq : data){
			SequenceAnnotation[] anns = seq.getAnnotation();
			for(int i=0;anns != null && i<anns.length;i++){
				if(anns[i] instanceof MotifAnnotation && (motifName == null || motifName.equals( anns[i].getIdentifier() ) ) ){
					if( ( ((MotifAnnotation)anns[i]).getStrandedness() == Strand.REVERSE && !rc) ||
							( ((MotifAnnotation)anns[i]).getStrandedness() != Strand.REVERSE && rc) ){
						seqs.add( seq.reverseComplement( ((MotifAnnotation)anns[i]).getPosition(), ((MotifAnnotation)anns[i]).getEnd() ) );
					}else{
						seqs.add( seq.getSubSequence( ((MotifAnnotation)anns[i]).getPosition(), ((MotifAnnotation)anns[i]).getLength() ) );
					}
				}
			}
		}
		return new DataSet("excised sites",seqs.toArray( new Sequence[0] ));
	}
	
	public static void main(String[] args) throws Exception {

		String[] prefixes = new String[]{"MA0001.1","MA0005.1","MA0015","MA0048","MA0048_ma52gauss","MA0048_ma52unif","MA0052","MA0054","MA0077"};//,"MA0001.1_in_hs"};
		String[] dataSet = new String[]{"MA0001","MA0005","MA0015","MA0048","MA0048+MA0052(gauss)","MA0048+MA0052(unif)","MA0052","MA0054","MA0077"};//,"MA0001.1_in_hs"};
		String[] orgs = new String[]{"at","at","dm","hs","hs","hs","hs","petunia","hs","at"};
		String[] infs = new String[]{"gauss","unif"};
		int[] lengths = new int[]{10,10,9,10,10,10,10,9,8,10};

		REnvironment r = null;
		try {
			r = new REnvironment( host, user, pass );
			
			double[][][][][] ppvs = new double[infs.length][prefixes.length][2][][];
			
			double[] sns = new double[]{0.1,0.25,0.3,0.5,0.7,0.75,0.9};
			String[] names = null;
			
			PFMDistance[] dists = new PFMDistance[]{
					//new PFMComparator.NormalizedEuclideanDistance(),
					//new PFMComparator.MinusPearsonCorrelationCoefficientPlusOne(),
					new PFMComparator.SymmetricKullbackLeiblerDivergence( 1 )
				};
			
			for(int pr=0;pr<prefixes.length;pr++){
				for(int in=0;in<infs.length;in++){
					for(int gl=0;gl<2;gl++){
						Discoverer.resetUsed();
						AlphabetContainer con = DNAAlphabetContainer.SINGLETON;
						String prefix = prefixes[pr];
						String org = orgs[pr];
						String inf = infs[in], infix = inf + "";
						int length = lengths[pr];
						boolean givenLength = gl==1;
						System.out.println("##################################################################################");
						System.out.println("data/" + DATATYPE +prefix+".sites_"+org+".txt_e.txt_"+inf+".txt");
						DataSet data = new DataSet( con, new SparseStringExtractor( HOME + "data/artificial/"+prefix+".sites_"+org+".txt_e.txt_"+inf+".txt" ) );
						System.out.println(HOME + "/data/" + DATATYPE + "/pos_"+prefix+".sites_"+org+".txt_e.txt_"+inf+".txt");
						DataSet truth = addPosition( data, HOME + "/data/" + DATATYPE + "/pos_"+prefix+".sites_"+org+".txt_e.txt_"+inf+".txt_new.txt",
								HOME + "/data/" + DATATYPE + "/revcom_"+prefix+".sites_"+org+".txt_e.txt.txt", length );
						
						double[] bg = PFMComparator.getCounts( truth );
						double[][] truePWMfw = getPWM( truth, motifName, false );
						double[][] truePWMrev = getPWM( truth, motifName, true );
						plotSeqLogo( r, HOME+"results/" + DATATYPE + "seqLogos/seqLogo-truePWM-"+prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+".pdf", truePWMfw );
	
						//aggregate results of other tools:
						LinkedList<DiscoveryResult> list = new LinkedList<DiscoveryResult>();
	
						//list.add( new DiscoveryResult( Discoverer.A__GLAM, truth, data, "mads2" + infix + "_anchor500_10.txt" ) );
						if(!givenLength){
							list.add( new DiscoveryResult( Discoverer.A__GLAM, truth, data, true, prefix+"_"+infix+"_anchor500.txt" ) );
							//list.add( new DiscoveryResult( Discoverer.A__GLAM, truth, data, prefix+"_"+infix+"_noanchor.txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.A__GLAM, truth, data, true, prefix+"_"+infix+"_anchor500_w"+length+".txt" ) );
							//list.add( new DiscoveryResult( Discoverer.A__GLAM, truth, data, prefix+"_"+infix+"_noanchor_w"+length+".txt" ) );
						}
	
						//	list.add( new DiscoveryResult( Discoverer.AMADEUS, truth, data, "targets-10-mads2" + infix + "-AT.txt", "" + 10 ) );
						//	list.add( new DiscoveryResult( Discoverer.AMADEUS, truth, data, "targets-15-mads2" + infix + "-AT.txt", "" + 15 ) );
	
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.DEME, truth, data, true, prefix+"_"+infix+"_w"+length+".txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.DEME, truth, data, true, prefix+"_"+infix+"_w15.txt" ) );
						}
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.DME, truth, data, true, "dme_"+prefix+"_"+infix+"_w"+length+".txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.DME, truth, data, true, "dme_"+prefix+"_"+infix+"_w15.txt" ) );
						}
	
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.GIBBS_SAMPLER, truth, data, true, prefix+"_"+infix+"_w"+length+".txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.GIBBS_SAMPLER, truth, data, true, prefix+"_"+infix+"_w15.txt" ) );
						}
						/*
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.CENTROID_GIBBS_SAMPLER, truth, data, true, prefix+"_"+infix+"_w"+length+"_centroid.txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.CENTROID_GIBBS_SAMPLER, truth, data, true, prefix+"_"+infix+"_w15_centroid.txt" ) );
						}/**/
							
						if(!givenLength){
							list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, true, prefix+"_"+infix+"_improbizer.txt", prefix+"_"+infix+"_motif.txt" ) );
							//list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, prefix+"_"+infix+"_improbizer_mo.txt", prefix+"_"+infix+"_motif.txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, true, prefix+"_"+infix+"_improbizer.txt", prefix+"_"+infix+"_motif_w"+length+".txt" ) );
							//list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, prefix+"_"+infix+"_improbizer_w"+length+"_mo.txt", prefix+"_"+infix+"_motif_w"+length+".txt" ) );
						}
						//	list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, "hits_mads2" + infix + "_w10_maxOcc3.txt", "improb_mads2" + infix + "_10.txt" ) );
						//	list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, "hits_mads2" + infix + "_maxOcc1.txt", "improb_mads2" + infix + ".txt" ) );
						//	list.add( new DiscoveryResult( Discoverer.IMPROBIZER, truth, data, "hits_mads2" + infix + "_maxOcc3.txt", "improb_mads2" + infix + ".txt" ) );
	
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.MEME, truth, data, true, prefix+"_"+infix+"_w"+length+".txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.MEME, truth, data, true, prefix+"_"+infix+"_w6-20.txt" ) );
						}
						
						/*TODO
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.MOAN, truth, data, true, "moan_fl_.*_"+prefix+"\\.sites_"+org+"\\.txt_e\\.txt_"+infix+"\\.txt" ) );
						}else{
							list.add( new DiscoveryResult( Discoverer.MOAN, truth, data, true, "moan_vl_.*_"+prefix+"\\.sites_"+org+"\\.txt_e\\.txt_"+infix+"\\.txt" ) );
						}*/
						
						if(givenLength){
							list.add( new DiscoveryResult( Discoverer.WEEDER, truth, data, true, "fl/" + prefix+".sites_"+org+".txt_e.txt_"+infix+".txt.fasta_fl.wee" ) );
						} else {
							list.add( new DiscoveryResult( Discoverer.WEEDER, truth, data, true, "vl/" + prefix+".sites_"+org+".txt_e.txt_"+infix+".txt.fasta.wee" ) );
						}
						//ParserToolBox.annotateWithWeederResults( data, false, home + "weeder/mads2" + infix + "_best.txt", ), //TODO
						
						//Sample bgSample = new Sample( con, new SparseStringExtractor( HOME + "data/" + DATATYPE + "/"+prefix+".sites_"+org+".txt_e.txt_bgonly.txt" ) );

						//list.add( new DiscoveryResult( Discoverer.DiPoMM, truth, data, false, getFiles( givenLength, prefix, infix, length, 75 )[0].file.getAbsolutePath(), HOME + "data/" + DATATYPE + "/"+prefix+".sites_"+org+".txt_e.txt_bgonly.txt", "1E-4" ) );
						//list.add( new DiscoveryResult( Discoverer.DiPoMM, truth, data, false, getFiles( givenLength, prefix, infix, length, 100 )[0].file.getAbsolutePath(), HOME + "data/" + DATATYPE + "/"+prefix+".sites_"+org+".txt_e.txt_bgonly.txt", "1E-4" ) );
						//list.add( new DiscoveryResult( Discoverer.DiPoMM, truth, data, false, getFiles( givenLength, prefix, infix, length, 125 )[0].file.getAbsolutePath(), HOME + "data/" + DATATYPE + "/"+prefix+".sites_"+org+".txt_e.txt_bgonly.txt", "1E-4" ) );
						FileAndName[] fnn = getFiles( givenLength, prefix, infix, length, "dipomm-150" );
						if( fnn != null && fnn.length > 0 ) {
							list.add( new DiscoveryResult( Discoverer.Dispom, truth, data, false, fnn[0].file.getAbsolutePath(), HOME + "data/" + DATATYPE + "/"+prefix+".sites_"+org+".txt_e.txt_bgonly.txt", "1E-4" ) );
						}
						/*
						fnn = getFiles( givenLength, prefix, infix, length, "dispom" );
						if( fnn != null && fnn.length > 0 ) {
							list.add( new DiscoveryResult( Discoverer.Dispom, truth, data, false, fnn[0].file.getAbsolutePath(), HOME + "data/" + DATATYPE + "/"+prefix+".sites_"+org+".txt_e.txt_bgonly.txt", "1E-4" ) );
						}*/
						
						Iterator<DiscoveryResult> it = null;
						/*PrintWriter wr = new PrintWriter(HOME+"/measures-"+prefix.replaceAll( "[\\._\\s]", "-" )+"_"+infix+"_givenLength-"+givenLength+".txt");
						wr.println( list.get(0).getHeading() );
						Iterator<DiscoveryResult> it = list.iterator();
						while( it.hasNext() ) {
							wr.println( it.next() );
						}
						wr.println();*/
						
						ppvs[in][pr][gl] = new double[list.size()][sns.length];
						names = new String[list.size()];
						it = list.iterator();
						int p = 0, fw = 0;
						ComparableElement<Boolean, Double> pwmComparison;
						PrintWriter wr2;
						while(it.hasNext()){
							DiscoveryResult dr = it.next();
							for(int s=0;s<sns.length;s++){
								ppvs[in][pr][gl][p][s] = dr.getPPVFor( sns[s] );
							}
							
							double[][] pwm = dr.getPWM();
							
							wr2 = new PrintWriter( HOME+"results/" + DATATYPE + "distances/dist-"+dr.getName().replaceAll( "[\\._\\s]", "-" )+"-"+prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+"-givenLength-"+givenLength+".txt" );
							fw = 0;
							for(int d=0;d<dists.length;d++){
								pwmComparison = getMeanDivergence( bg, pwm, truePWMfw, truePWMrev, dists[d] );
								wr2.println(/*dists[d].getClass().getSimpleName()+":\t"*/+ pwmComparison.getWeight() /*+"\t\\\\"*/ );
								if( pwmComparison.getElement() ) {
									fw++;
								}
							}
							wr2.close();
							
							dr.createPlots( r, fw >= dists.length/2d, prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+"-givenLength-"+givenLength );
							//plotSeqLogo( r, HOME+"results/" + DATATYPE + "seqLogos/seqLogo-"+dr.getName().replaceAll( "[\\._\\s]", "-" )+"-"+prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+"-givenLength-"+givenLength+".pdf", fw >= dists.length/2d? pwm : PFMComparator.getReverseComplement((DNAAlphabet)con.getAlphabetAt(0), pwm) );
							//System.out.print(dr.getName()+"\t");
							
							names[p] = dr.getName();
							p++;
						}					
						
						//remove bad results
						DiscoveryResult dr;
						for( int i=0; i < list.size(); i++ ) {
							dr = list.get(i);
							if( !dr.plot() ) {
								list.remove(i--);
							}
						}	
						
						LinkedList<DoubleTableResult> dtrPR = new LinkedList<DoubleTableResult>();
						LinkedList<DoubleTableResult> dtrROC = new LinkedList<DoubleTableResult>();
						ArrayList<String> colors = new ArrayList<String>();
							
						StringBuffer points = new StringBuffer( list.size()*100 );
						points.append( "\n" );
						String n;
						for( int i = 0; i < list.size(); i++ ) {
							dr = list.get( i );
							points.append( "points( " + dr.getValue(0) + ", " + dr.getValue(2) + ", col=\"" + dr.getColor()
									+ "\", pch=" +dr.getPCH() + ", cex=2, lwd=3 );\n"  );
							n = dr.getName();
							points.append( "text( " + dr.getValue(0) + ", " + dr.getValue(2) + ", \"" + n + "\", col=\"" + dr.getColor() + "\", adj=c(1,1.5), cex=1.75 );\n"  );
							dtrPR.add( dr.getPrCurve() );
							dtrROC.add( dr.getROCCurve() );
							colors.add( dr.getColor() );
						}
						//System.out.println( "plot(0,0,xlim=c(0,1),ylim=c(0,1),col=0)" );
						//System.out.println( points );
	
						/*
						//DiPoMM results
						names[p] = "DiPoMM";
						ScoreClassifier cl;
						MotifDiscoverer md;
	
						double sign = 1E-3, auc;
						SignificantMotifOccurrencesFinder smof;
	
						LinkedList<double[]> listOfDoublesPR = new LinkedList<double[]>();
						LinkedList<double[]> listOfDoublesROC = new LinkedList<double[]>();
						double[][] scores;
						double[][] sortedScores;
	
						//cl = new GenDisMixClassifier( FileManager.readFile( new File( HOME + "/results/our_tool/classifier-single-mads2" + infix + ".txt-bg" + bg + ".txt-order=0-0-init=15.xml" ) ) );
						FileAndName[] files = getFiles( givenLength, prefix, infix, length );
						int c=1;
						for( FileAndName f: files ) {
							System.out.println(HOME + "data/" + DATATYPE +prefix+".sites_"+org+".txt_e.txt_"+inf+".txt -> "+f.file.getAbsolutePath());
							
							//wr.println( f.getName() );
							try{
								cl = new MSPClassifier( FileManager.readFile( f.file ) );
							}catch(NonParsableException e){
								cl = new GenDisMixClassifier( FileManager.readFile( f.file ) );
							}
							md = (MotifDiscoverer) cl.getScoringFunction( 0 );
	
							smof = new SignificantMotifOccurrencesFinder(md,bgSample,sign);
							//	smof = new SignificantMotifOccurrencesFinder(md,RandomSeqType.PERMUTED,1000,sign);
							for( int i = 0; i < 1; i++ ) {
								listOfDoublesPR.clear();
								scores = smof.getValuesForEachNucleotide( data, 0, 0, i==1 );
								sortedScores = MotifDiscoveryAssessment.getSortedValuesForMotifAndFlanking( truth, scores, smof.getOffsetForAucPR(), smof.getFactorForAucPR(), motifName );
								auc = ScoreBasedPerformanceMeasureDefinitions.getAUC_PR( sortedScores[0], sortedScores[1], listOfDoublesPR );
								ScoreBasedPerformanceMeasureDefinitions.getAUC_ROC( sortedScores[0], sortedScores[1], listOfDoublesROC );
								
								Sample pred = smof.annotateMotifs( data );
								double[][] pwm = getPWM( pred, null, false );
								wr2 = new PrintWriter( HOME+"results/" + DATATYPE + "distances/dist-"+f.name.replaceAll( "[\\._\\s]", "-" )+"-"+prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+"-givenLength-"+givenLength+".txt" );
								
								fw=0;
								for(int d=0;d<dists.length;d++){
									pwmComparison = getMeanDivergence( bg, pwm, truePWMfw, truePWMrev, dists[d] );
									wr2.println(//dists[d].getClass().getSimpleName()+":\t"+
										pwmComparison.getWeight()//+"\t\\\\"
									);
									if( pwmComparison.getElement() ) {
										fw++;
									}
								}
								wr2.close();
								
								plotSeqLogo( r, HOME+"results/" + DATATYPE + "seqLogos/seqLogo-"+f.name.replaceAll( "[\\._\\s]", "-" )+"-"+prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+"-givenLength-"+givenLength+".pdf", fw >= dists.length/2d? pwm : PFMComparator.getReverseComplement((DNAAlphabet)con.getAlphabetAt(0), pwm) );
								
								plotPositionalDistribution( r, HOME+"results/" + DATATYPE + "positionalDistribution/pos-"+f.name.replaceAll( "[\\._\\s]", "-" )+"-"+prefix.replaceAll( "[\\._\\s]", "-" )+"-"+infix+"-givenLength-"+givenLength+".pdf", (ExtendedZOOPSDiffSM) md, truth );
								
								System.out.println("DiPoMM: " + auc);
								dtrPR.add( new DoubleTableResult( Measure.PrecisionRecallCurve.getNameString(), Measure.PrecisionRecallCurve.getCommentString(), listOfDoublesPR ) );
								dtrROC.add( new DoubleTableResult( Measure.ReceiverOperatingCharacteristicCurve.getNameString(), Measure.ReceiverOperatingCharacteristicCurve.getCommentString(), listOfDoublesROC ) );
								points.append( "text( " + 1.0 + ", " + (1.0-(c*0.05)) + ", \"" + (f.name) + "\", col=" + c + ", pos=2);\n"  );
								colors.add( c );
	
								for(int s=0;s<sns.length;s++){
									ppvs[in][pr][gl][ ppvs[in][pr][gl].length-1 ][s] = getPPVFor( dtrPR.getLast(), sns[s] );
								}	
								c++;
							}
						}
						/**/
						//wr.close();

						//System.out.println("DTR: "+dtr.size()+" COLORS: "+colors.length());
						System.out.println( colors );
						StringBuffer plotcmd =  DoubleTableResult.getPlotCommands( r, ", xlim=c(0, 1), ylim=c(0, 1), xlab=\"Nucleotide Recall\", ylab=\"Nucleotide Precision\", lwd=3", colors.toArray(new String[0]),dtrPR.toArray( new DoubleTableResult[0] ) );
						//System.out.println(plotcmd);
						String palette = "par(mar=c(5,5,0,0),cex.lab=2,cex.main=2,cex.axis=2);\n";
						//REnvironment.showImage( "pr-"+prefix.replaceAll( "[\\._]", "-" )+"_"+infix+"_givenLength-"+givenLength, r.plot( palette + plotcmd.toString() + "\n" + points.toString() ), JFrame.EXIT_ON_CLOSE );
						r.plotToPDF( palette+plotcmd.toString()+ "\n" + points.toString(), HOME+"results/" + DATATYPE + "pics/plot-pr-"+prefix.replaceAll( "[\\._]", "-" )+"_"+infix+"_givenLength-"+givenLength+".pdf", true );
						
						//plotcmd = DoubleTableResult.getPlotCommands( r, null, colors.toArray(), dtrROC.toArray( new DoubleTableResult[0] ) );
						//r.plotToPDF( palette+plotcmd.toString() /*+ "\n" + points.toString()*/, HOME+"results/" + DATATYPE + "pics/plot-roc-"+prefix.replaceAll( "[\\._]", "-" )+"_"+infix+"_givenLength-"+givenLength+".pdf", true );
						/**/
					}
				}
			}
			
			
			for(int i=0;i<infs.length;i++){
				for(int j=0;j<2;j++){
					for(int k=0;k<sns.length;k++){
						PrintWriter pr = new PrintWriter( HOME+"results/" + DATATYPE + "tables/table_"+infs[i]+"_"+(j==0 ? "variable" : "fixed")+"_"+sns[k]+".txt" );
						
						for(int l=0;l<prefixes.length;l++){
							pr.print( "\t\""+dataSet[l] + "\"" );
						}
						pr.println();						
						for(int l=0;l<ppvs[i][0][j].length;l++){
							pr.print( names[l] );
							for(int m=0;m<prefixes.length;m++){
								if(Double.isNaN( ppvs[i][m][j][l][k] )){
									pr.print( "\t0" );//TODO pr.print( "\tNA" );
								}else{
									pr.print( "\t"+ppvs[i][m][j][l][k] );
								}
							}
							pr.println();
						}
						
						pr.close();
					}
				}
			}
		}catch( Exception e ) {
			e.printStackTrace();
		} finally {
			if( r != null ) {
				r.close();
				System.out.println( "closed REnvironment" );
			}
		}		
	}
	
	
	private static FileAndName[] getFiles(boolean givenLength, String prefix, String infix, int length, String tool ){
		LinkedList<File> files = new LinkedList<File>();
		LinkedList<String> names = new LinkedList<String>();
		/*File[] temp = (givenLength ? new File( HOME+"results/" + DATATYPE + TOOLS+"/dipomm/fixedlength/" ).listFiles( (FileFilter)
				new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-length-"+length+"-adjust-false.*" ) )
		                             :new File( HOME+"results/" + DATATYPE + TOOLS+"/dipomm/variablelength/" ).listFiles( (FileFilter)
				new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-adjust-true.*" ) ) );
		
		for(int i=0;i<temp.length;i++){
			files.add( temp[i] );
			names.add( "DiPoMM" );
		}
		temp = (givenLength ? new File( HOME+"results/" + DATATYPE + TOOLS+"/dipomm/ord_0-cut/" ).listFiles( (FileFilter)
				new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-length-"+length+"-adjust-false.*" ) )
		                             :new File( HOME+"results/" + DATATYPE + TOOLS+"/dipomm/ord_0-cut/" ).listFiles( (FileFilter)
				new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-adjust-true.*" ) ) );
		
		for(int i=0;i<temp.length;i++){
			files.add( temp[i] );
			names.add( "DiPoMM cut" );
		}
		temp = (givenLength ? new File( HOME+"results/" + DATATYPE + TOOLS+"/dipomm/ord_0-shift/" ).listFiles( (FileFilter)
				new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-length-"+length+"-adjust-false.*" ) )
		                             :new File( HOME+"results/" + DATATYPE + TOOLS+"/dipomm/ord_0-shift/" ).listFiles( (FileFilter)
				new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-adjust-true.*" ) ) );
		
		for(int i=0;i<temp.length;i++){
			files.add( temp[i] );
			names.add( "DiPoMM shift" );
		}*/
		
		FileFilter regExFF = givenLength
		    ? new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-length-"+length+"-adjust-false.*" ) 
			: new RegExFilenameFilter("",false,true,".*classifier-"+prefix+"\\..*_"+ infix +"\\..*-adjust-true.*" );
		File[] temp = new File( HOME+"results/" + DATATYPE + TOOLS+"/" + tool + "/" ).listFiles( regExFF );
		
		if( temp.length == 0 ) {
			System.out.println( regExFF );
		}
		
		for(int i=0;i<temp.length;i++){
			files.add( temp[i] );
			names.add( tool );//"Dispom" );
		}
		
		FileAndName[] fnn = new FileAndName[files.size()];
		for(int i=0;i<fnn.length;i++){
			fnn[i] = new FileAndName(files.get( i ),names.get( i ));
		}
		return fnn;
	}
	
	private static class FileAndName{
		
		public File file;
		public String name;
		
		public FileAndName(File file, String name){
			this.file = file;
			this.name = name;
		}
		
	}
}
