package projects.dream2016;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.broad.igv.bbfile.BBFileHeader;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import de.jstacs.utils.Time;
import de.jstacs.utils.ToolBox;


public class BigWig2Interval {
	
	/**
	 * Creates genomewide DNase interval files.
	 * 
	 * @param args
	 * 0 .. bigwig input
	 * 1 .. reference file
	 * 2 ... blacklistfiltered
	 * 3 ... boolean: true (=simple) or false (=orange,...)
	 * 
	 * @throws IOException
	 */
	public static void main( String[] args ) throws IOException {
		BigWigAccessor bwa = new BigWigAccessor(args[0]);
		GZIPInputStream stream = new GZIPInputStream(new FileInputStream(args[1]));
		BufferedReader read = new BufferedReader(new InputStreamReader(stream));
		boolean simple = Boolean.parseBoolean(args[3]);
		int anz = simple ? 4 : 5;
		
		//blacklistfiltered
		File f = new File( args[2] );
		BufferedReader r = new BufferedReader( new FileReader(f) );
		HashMap<String, ArrayList<int[]>> hash = new HashMap<String, ArrayList<int[]>>();
		String line;
		while( (line = r.readLine()) != null ) {
			String[] split = line.split("\t");
			ArrayList<int[]> current = hash.get(split[0]);
			if( current == null ) {
				current = new ArrayList<int[]>();
				hash.put(split[0], current);
			}
			current.add(new int[]{Integer.parseInt(split[1]), Integer.parseInt(split[2])} );
		}
		r.close();
		
		String chr=null;
		int s = -100, e;
		//double[] profile = new double[50];
		
		FileOutputStream output = new FileOutputStream(args[0]+"-interval" + (simple?"":"-orange") + ".txt.gz");
		try {
			BufferedWriter w = new BufferedWriter( new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8") );
			Time t = Time.getTimeInstance(null);
			try {
				/*
				while( (line=read.readLine()) != null ) {
					if(line.startsWith("[")){
						System.out.println(line);
						chr = line.substring(1, line.length()-1);
						s=0;
						w.append(line);
						w.newLine();
					}else{
						bwa.fillProfileInRegion(chr, s,s+50,profile);
						w.append(ToolBox.min(profile) + "\t" + ToolBox.median(profile) + "\t" + ToolBox.max(profile));
						w.newLine();
						s+=50;
					}
				}
				*/
				
				line = read.readLine();
				do {
					System.out.println(line);
					chr = line.substring(1, line.length()-1);
					ArrayList<int[]> interval = hash.get(chr);
					int i = 0;
					int[] inter=null;
					s=e=0;
					if( !simple ) {
						w.append(line); w.newLine();
					}
					
					while( (line=read.readLine()) != null && line.charAt(0)!='[' ) {
						e+=50;
					}
					double[] profile = bwa.getProfileInRegion(chr, s, e);
					/*
					System.out.println(chr);
					for( int i = 1000000; i <1100000; i++) {
						System.out.println(profile[i]);
					}/**/
					while( s+50<=e ) {
						if( simple ){
							w.append(chr + "\t" + s + "\t" );//TODO
						}
						while( interval != null && i < interval.size() && s >= (inter=interval.get(i))[1] ) {
							i++;
						}
						if( interval != null && i < interval.size() && s >= inter[0] ) {
							if( simple ) {
								w.append( ToolBox.min(s,s+50,profile) + "\t" + ToolBox.median(s,s+50,profile) + "\t" + ToolBox.max(s,s+50,profile) + "\t" + ToolBox.percentile(s, s+50, profile, 0.25) );  
							} else {
								//TODO interval size: autosome 400bp -> JTeam 450bp?
								w.append( orange( s-200, s+250, profile )
									+ "\t" + mostMonotonSteps( s-200, s, profile, 1 )
									+ "\t" + mostMonotonSteps( s-200, s, profile, -1 )
									+ "\t" + mostMonotonSteps( s+50,s+250, profile, 1 )
									+ "\t" + mostMonotonSteps( s+50, s+250, profile, -1 ) );
							}
						} else {
							for( int j = 0; j < anz; j++ ) {//TODO length
								w.append((j==0?"":"\t") + 0);
							}
						}
						w.newLine();
						s=s+50;
					}					
					System.out.println(t.getElapsedTime());
				} while( line != null );
			} finally {
				w.close();
			}
		} finally {
			output.close();
		}
		read.close();
		stream.close();
	}
	
	private static int orange( int start, int end, double[] profile ) {
		int different = 0;
		start = Math.max(0, start);
		end = Math.min(profile.length, end);
		for( int s = start+1; s < end; s++ ) {
			if( profile[s-1] != profile[s] ) {
				different++;
			}
		}
		return different;
	}
	
	private static int mostMonotonSteps( int start, int end, double[] profile, double vz ) {
		int num=0, max=0;
		start = Math.max(0, start);
		end = Math.min(profile.length, end);
		for( int s = start+1; s < end; s++ ) {
			if( vz*profile[s-1] < vz*profile[s] ) {
				num++;
			} else if( vz*profile[s-1] > vz*profile[s] ) {
				if( num > max ) {
					max=num;
				}
				num=0;
			}
		}
		return max;
	}	
	
	private static class BigWigAccessor {

		private BBFileReader reader;
		
		public BigWigAccessor(String bigWigFile) throws IOException {
			reader = new BBFileReader(bigWigFile);
			
			BBFileHeader header = reader.getBBFileHeader();
			
			if(!header.isHeaderOK()){
				throw new RuntimeException("Header not OK");
			}
			
			if(!header.isBigWig()){
				throw new RuntimeException("No Bigwig");
			}
		}
		
		
		public double[] getProfileInRegion(String chr, int start, int end){
			double[] res = new double[end-start];
			fillProfileInRegion( chr, start, end, res );
			return res;
		}
		
		public void fillProfileInRegion(String chr, int start, int end, double[] res) {
			BigWigIterator it = reader.getBigWigIterator(chr,start,chr,end,false);
			
			while(it.hasNext()){
				WigItem item = it.next();
				Arrays.fill( res, Math.max(0, item.getStartBase()-start), Math.min(item.getEndBase(),end)-start, item.getWigValue() );
				/*
				int s = item.getStartBase();
				int e = item.getEndBase();
				double v = item.getWigValue();
				for(int i=Math.max(start, s);i<Math.min(e, end);i++){
					//System.out.println(i+"\t"+(i-start)+"\t"+v);
					res[i-start] = v;
				}*/
			}		
		}	
	}
	
}


