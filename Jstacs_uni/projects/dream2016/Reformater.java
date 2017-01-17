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

public class Reformater {

	public static void main(String[] args) throws IOException {
		BufferedReader r = new BufferedReader( new InputStreamReader(new GZIPInputStream(new FileInputStream(args[0]))) );
		BufferedWriter w = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(args[0]+"-r.gz"))) );
		String chr = null, line;
		String[] split;
		while( (line=r.readLine()) != null ) {
			split = line.split("\t");
			if( chr == null || !split[0].equals(chr) ) {
				chr = split[0];
				w.append("["+chr+"]");
				w.newLine();
			}
			for( int j = 2; j < split.length; j++ ) {
				w.append( (j==2?"":"\t") + split[j] );
			}
			w.newLine();
		}
		r.close();
		w.close();
	}
}