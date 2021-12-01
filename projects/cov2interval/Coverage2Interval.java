import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 * Simple class for parsing coverage files.
 * 
 * @author Jens Keilwagen
 */
public class Coverage2Interval {

	/**
	 * A class that stores the new chromosome name and the offset
	 * 
	 * @author Jens Keilwagen
	 */
	static class Remapping {
		String chr;
		int offset;
		
		Remapping( String newChr, String o ) {
			chr = newChr;
			offset=Integer.parseInt(o);
		}
	}
	
	public static void main(String[] args) throws Exception {
		if( args.length==0 ) {
			System.out.println("command line: java Coverage2Interval <read mapping> <chromosome mapping> [<maximal coverage>]");
			System.exit(0);
		}
		
		BufferedReader r;
		String line;
		
		//read chromosome remapping
		HashMap<String,Remapping> hash = null;
		File f = new File(args[1]);
		if( f.exists() ) {
			hash = new HashMap<String,Remapping>();
			r = new BufferedReader( new FileReader(args[1]) );
			while( (line=r.readLine()) != null ) {
				String[] split = line.split("\t");
				hash.put(split[0], new Remapping(split[1],split[2]));
			}
			r.close();
		}
		
		//open read mapping
		if( args[0].equals("-") ) {
			r = new BufferedReader( new InputStreamReader(System.in) );
		} else {
			r = new BufferedReader( new FileReader(args[0]) );
		}
		String oldChr=null;
		long startPos=-100, lastPos=-99;
		int oldCov=-100;
		int maxCov = args.length==2 ? Integer.MAX_VALUE : Integer.parseInt(args[2]);
		if( maxCov<=0 ) throw new IllegalArgumentException("The maximal covergae has to be positive.");
		Remapping rm = null;
		while( (line=r.readLine()) != null ) {
			String[] split = line.split("\t");
			String newChr=split[0];
			int pos = Integer.parseInt(split[1]);
			if( hash != null && (rm=hash.get(newChr))!=null ) {
				newChr=rm.chr;
				pos +=rm.offset;
			}
			//sum coverage (if necessary)
			int cov = 0;
			for( int i=2; i < split.length; i++ ) {
				cov += Integer.parseInt(split[i]);
			}			
			//replace with maxCoverage (if necessary)
			if( cov >= maxCov ) {
				cov=maxCov;
			}
System.out.println("TBR\t" + line + "\t" + newChr + "\t" + pos + "\t" + cov);

			//output
			if( !newChr.equals(oldChr) ) {
				//new chromosome
				if( oldChr!=null ) System.out.println( oldChr + "\t" + startPos + "\t" + lastPos + "\t" + oldCov );
				oldChr = newChr;
				startPos = lastPos = pos;
				oldCov = cov;
			} else {
				if( lastPos+1 == pos ) {
					//adjacent positions
					if( oldCov==cov ) {
						lastPos=pos;
					} else {
						System.out.println( oldChr + "\t" + startPos + "\t" + lastPos + "\t" + oldCov );
						oldChr = split[0];
						startPos = lastPos=pos;
						oldCov = cov;
					}
				} else {
					//new region
					System.out.println( oldChr + "\t" + startPos + "\t" + lastPos + "\t" + oldCov );
					//System.out.println( oldChr + "\t" + (lastPos+1) + "\t" + (pos-1) + "\t" + 0 );
					oldChr = split[0];
					startPos = lastPos=pos;
					oldCov = cov;
				}
			}
		}
		//last region
		System.out.println( oldChr + "\t" + startPos + "\t" + lastPos + "\t" + oldCov );
		r.close();
	}
}
