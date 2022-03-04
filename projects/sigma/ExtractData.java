package projects.sigma;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;

import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.sequences.Sequence;

public class ExtractData {

	public static void main(String[] args) throws Exception {
		
		String coords = args[0];
		HashMap<String,DataSet> chrs = new HashMap<String,DataSet>();
		for(int i=1;i<args.length;i++){
			String key = args[i].replaceAll("\\.fasta", "").replaceAll(".*/", "");
			chrs.put(key, new DNADataSet(args[i]) );
		}

		BufferedReader read = new BufferedReader(new FileReader(coords));
		
		String str = null;
		
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			String chr = parts[1].replaceAll("[\\-\\+]$", "");
			int pos = Integer.parseInt(parts[2]);
			boolean fwd = parts[3].equals("+");
			
			Sequence seq = chrs.get(chr).getElementAt(0);
			Sequence sub = null;
			if(fwd){
				if(pos-90<0 || pos+10>seq.getLength()) {
					System.err.println("EX: "+parts[0]+" "+chr);
					continue;
				}
				sub = seq.getSubSequence(pos-90, 100);
			}else{
				if(pos-10<0 || pos+90>seq.getLength()) {
					System.err.println("EX: "+parts[0]+" "+chr);
					continue;
				}
				sub = seq.getSubSequence(pos-10, 100).reverseComplement();
			}
			
			System.out.println("> "+parts[0]+" "+chr+" "+pos+" "+fwd+"\n"+sub);
			
		}
		read.close();
		
	}

}
