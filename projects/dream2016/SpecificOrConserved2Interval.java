package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import cern.colt.Arrays;
import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.RegExFilenameFilter.Directory;
import de.jstacs.utils.ToolBox;

/**
 * 
 * @author Jens Keilwagen
 */
public class SpecificOrConserved2Interval {

	/**
	 * 
	 * @param args
	 * 0 reg. expression
	 * 1 home
	 * 2 columns (comma separated)
	 * 3 conserved prefix
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		FilenameFilter filter = new RegExFilenameFilter("", Directory.FORBIDDEN, true, args[0] );
		String[] fName = new File( args[1] ).list(filter);
		String[] s = args[2].split(",");
		int[] cols = new int[s.length];
		for( int i = 0; i < cols.length; i++ ) {
			cols[i] = Integer.parseInt( s[i] );
		}
		System.out.println(Arrays.toString(cols));
		
		BufferedReader[] reader = new BufferedReader[fName.length];
		BufferedWriter[] writer = new BufferedWriter[fName.length+1];
		for(int f=0;f<reader.length;f++){
			System.out.println(fName[f]);
			
			GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[1]+"/"+fName[f]));
			reader[f] = new BufferedReader(new InputStreamReader(stream));
			
			FileOutputStream output = new FileOutputStream(args[1]+"/"+fName[f] + "-specific.txt.gz");
			writer[f] = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
		}
		FileOutputStream output = new FileOutputStream(args[1]+"/"+args[3] + "-conserved.txt.gz");
		writer[writer.length-1] = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
		double[] value = new double[reader.length];
		
		String[] line = new String[reader.length];
		String[][] split = new String[reader.length][];
		int a = 0;
		outerloop: while( true ) {
			//read
			boolean any = false;
			for(int f=0;f<reader.length;f++){
				line[f] = reader[f].readLine();
				if( line[f] == null ) {
					break outerloop;
				}
				any |= line[f].charAt(0) == '[';
			}
			//synchronize
			if( any ) {
				String chr = null;
				for(int f=0;f<reader.length;f++){
					if( line[f].charAt(0) != '[' ) {
						throw new RuntimeException( "File " + f + ": " + fName[f] );
					} else {
						while( line[f].charAt(0)!='[' || line[f].contains("_") ) {
							line[f] = reader[f].readLine();
						}
						if( line[f] != null ) {
							if( chr == null ) {
								chr = line[f];
							} else {
								if( !chr.equals( line[f] ) ) {
									throw new RuntimeException( "File " + f + ": " + fName[f] );
								}
							}
						} else {
							break outerloop;
						}
					}
				}
				System.out.println(chr);
				a=0;
			}
			
			//write
			if( line[0].charAt(0) == '[' ) {
				for(int f=0;f<writer.length;f++){
					writer[f].append(line[0]);
				}
			} else {
				//compute
				for(int f=0;f<reader.length;f++){
					split[f] = line[f].split("\t");
				}
				
				for( int c = 0; c < cols.length; c++ ) {
					for(int f=0;f<reader.length;f++){
						value[f] = Double.parseDouble(split[f][cols[c]]);
					}
					
					//specific
					double mean = ToolBox.mean(0,value.length,value);
					for(int f=0;f<reader.length;f++){
						writer[f].append( (c==0?"":"\t") + (mean==0?0:(value[f]/mean)) );
					}
					
					//conserved
					double sd = ToolBox.sd(0,value.length,value);
					writer[writer.length-1].append( (c==0?"":"\t") + ToolBox.min(value) + "\t" + (mean==0?0:(sd/mean)) );
				}
			}
			a++;
			for(int f=0;f<writer.length;f++){
				writer[f].newLine();
			}
			if( a%1000000 == 0 ) {
				System.out.println(a);
			}
		}
		
		for(int f=0;f<reader.length;f++){
			reader[f].close();
			writer[f].close();
		}
		writer[writer.length-1].close();
	}
}