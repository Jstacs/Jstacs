package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

import de.jstacs.classifiers.AbstractScoreBasedClassifier.DoubleTableResult;
import de.jstacs.classifiers.performanceMeasures.AucPR;
import de.jstacs.classifiers.performanceMeasures.AucROC;
import de.jstacs.classifiers.performanceMeasures.PRCurve;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.DoubleList;

public class Performance {

	/**
	 * @param args
	 * 0 .. agg file
	 * 1 .. labels
	 * 2 .. column of cell type
	 */
	public static void main(String[] args) throws IOException {
		boolean curve = false;
		
		GZIPInputStream stream2 = new GZIPInputStream(new FileInputStream(args[0]));
		BufferedReader preds = new BufferedReader(new InputStreamReader(stream2));
		
		int col=Integer.parseInt(args[2]);
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[1]));
		BufferedReader labels = new BufferedReader(new InputStreamReader(stream));
		
		String lab = labels.readLine(), pred;
		System.out.println( "Compute performance for cell type: " + lab.split("\t")[col] );
		
		DoubleList pos = new DoubleList();
		DoubleList neg = new DoubleList();
		
		String chr=null;
		while( (lab = labels.readLine()) != null && (pred = preds.readLine()) != null ) {
			if( chr == null || !lab.startsWith(chr) ) {
				chr = lab.substring(0,lab.indexOf('\t')+1 );
				System.out.println(chr);
			}
			/*System.out.println(pred);
			System.out.println(lab);
			System.exit(1);/**/
			//skip chr
			String old = null;
			while( pred != null && !pred.startsWith(chr) ) {
				if( old == null || !pred.startsWith(old) ) {
					old = pred.substring(0,pred.indexOf('\t')+1 );
					System.out.println("Skip " +old);
				}
				pred=preds.readLine();
			}
			
			if( pred == null ) {
				break;
			} else {
				String[] lSplit = lab.split("\t");
				String[] pSplit = pred.split("\t");
				if( lSplit[0].equalsIgnoreCase(pSplit[0]) && lSplit[1].equalsIgnoreCase(pSplit[1]) ) {
					double ag = Double.parseDouble(pSplit[3]); 
					switch( lSplit[col].charAt(0) ) {
						case 'b': case 'B': pos.add(ag);break;
						case 'u': case 'U': neg.add(ag);break;
					}
				} else {
					throw new RuntimeException("could not match:\n" + lab + "\n" + pred );
				}
			}
		}
		labels.close();
		preds.close();

		double[] po = pos.toArray();
		double[] ne = neg.toArray();
		Arrays.sort(po);
		Arrays.sort(ne);

		System.out.println( "#positives: " + po.length + "\t" + po[0]+ " .. " + po[po.length-1]);
		System.out.println( "#negatives: " + ne.length + "\t" + ne[0]+ " .. " + ne[ne.length-1] );
					
		System.out.println( "random: " + po.length / (double)(po.length + ne.length) );
		
		System.out.println();
		System.out.println( new AucROC().compute( po, ne ) );
		System.out.println( new AucPR().compute( po, ne ) );
		
		ResultSet rs = new PRCurve().compute(po, ne);
		DoubleTableResult dtr = (DoubleTableResult)rs.getResultAt(rs.getNumberOfResults()-1);
		double[][] c = dtr.getValue();
		double max10 = 0;
		double max50 = 0;
		BufferedWriter w = null;
		if( curve ) {
			w = new BufferedWriter(new FileWriter(args[0] +".prcurve"));
		}
		for( int i = 0; i < c.length; i++ ) {
			if( c[i][1] >= 0.9 && c[i][0] > max10 ) {
				max10 = c[i][0];
			}
			if( c[i][1] >= 0.5 && c[i][0] > max50 ) {
				max50 = c[i][0];
			}
			if( w != null ) {
				w.append( c[i][0] + "\t" + c[i][1] );
				w.newLine();
			}
		}
		if( w != null ) {
			w.close();
		}
		System.out.println( "Recall at 10% FDR: " + max10 );
		System.out.println( "Recall at 50% FDR: " + max50 );
	}
}