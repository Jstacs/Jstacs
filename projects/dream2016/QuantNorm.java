package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class QuantNorm {

	public static void main(String[] args) throws IOException {
		int col = Integer.parseInt(args[0]);
		BufferedReader[] r = new BufferedReader[args.length-1];
		BufferedWriter[] w = new BufferedWriter[args.length-1];
		for( int i = 0; i <r.length; i++ ) {
			r[i] = new BufferedReader( new InputStreamReader(new GZIPInputStream(new FileInputStream(args[1+i]))) );
			w[i] = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args[1+i]+"-q.gz"))) );
		}
		String[] line = new String[r.length];
		String[][] split = new String[r.length][];
		outerloop: do {
			double mean = 0;
			for( int i = 0; i <r.length; i++ ) {
				line[i] = r[i].readLine();
				if( line[i] != null && line[i].length() > 0 ) {
					split[i] = line[i].split("\t");
					mean += Double.parseDouble(split[i][col]);					
				} else {
					break outerloop;
				}
			}
			mean /= r.length;
			for( int i = 0; i <r.length; i++ ) {
				for( int j =0; j < split[i].length; j++ ) {
					w[i].append( (j==0?"":"\t") + (j==col?mean:split[i][j]) );
				}
				w[i].newLine();
			}
		} while( true );
		
		for( int i = 0; i <r.length; i++ ) {
			r[i].close();
			w[i].close();
		}
	}
}
