package projects;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

/**
 * This class splits fastA files in a user-specified number of fastA-files.
 * It helps to run tools in parallel for instance in Galaxy workflows.
 * 
 * @author Jens Keilwagen
 */
public class FastaSplitter {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		int num = Integer.parseInt(args[1]);
		BufferedReader r = new BufferedReader( new FileReader(args[0]) );
		String delim=null, id = null;
		HashMap<String, Integer> hash = null;
		if( args.length>2 ) {
			delim=args[2];
			hash = new HashMap<String, Integer>();
		}
		BufferedWriter[] w = new BufferedWriter[num];
		int[] stats = new int[num];
		for( int i = 0; i < num; i++ ) {
			w[i] = new BufferedWriter(new FileWriter("split-"+i+".fasta"));			
		}
		
		String line;
		int current = num-1;
		while( (line=r.readLine()) != null ) {
			
			if( line.length()>0 && line.charAt(0) == '>' ) {
				if( delim==null ) {
					current = (current+1) % num;
				} else {
					id = line.substring(1, line.lastIndexOf(delim) );
					System.out.println(id);
					Integer idx = hash.get(id);
					System.out.println(idx);
					if( idx == null ) {
						//this id occurs the first time => find a split with minimal number of sequences so far
						int min = 0;
						for( int i = 1; i < num; i++ ) {
							if( stats[i] < stats[min] ) {
								min = i;
							}
						}
						current = min;
						hash.put(id, current);
					} else {
						//this id was assigned to split idx before
						current = idx;
					}
				}				
				stats[current]++;
			}
			w[current].write(line);
			w[current].newLine();
		}
		r.close();
		for( int i = 0; i < num; i++ ) {
			w[i].close();	
			if( stats[i] == 0 ) {
				new File("split-"+i+".fasta").deleteOnExit();
			}
		}
	}

}
