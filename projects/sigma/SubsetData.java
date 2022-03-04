package projects.sigma;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;

import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;
import de.jstacs.utils.DoubleList;

public class SubsetData {

	public static void main(String[] args) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		
		DNADataSet data = new DNADataSet(args[0], '>', new SimpleSequenceAnnotationParser());
		
		HashMap<String,DoubleList> pos = new HashMap<>();
		
		BufferedReader read = new BufferedReader(new FileReader(args[1]));
		
		String str = read.readLine();
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			if(!pos.containsKey(parts[6])){
				pos.put(parts[6], new DoubleList());
			}
			pos.get(parts[6]).add(Double.parseDouble(parts[3]));
		}
		read.close();
		
		
		HashMap<String,Double> medians = new HashMap<>();
		for(String key : pos.keySet()){
			DoubleList temp = pos.get(key);
			medians.put(key, temp.median(0, temp.length()) );
		}
		
		
		read = new BufferedReader(new FileReader(args[1]));
		
		HashMap<String,LinkedList<Sequence>> map = new HashMap<>();
		
		str = read.readLine();
		while( (str = read.readLine()) != null ){
			String[] parts = str.split("\t");
			
			double l = Double.parseDouble(parts[3]);
			
			String key = parts[6];
			if(Math.abs(l-medians.get(key))<5){
				key += "p";
			}else{
				key += "n";
			}
			
			if(!map.containsKey(key)){
				map.put(key, new LinkedList<Sequence>());
			}
			map.get(key).add(data.getElementAt(Integer.parseInt(parts[0])-1));
			
		}
		read.close();
		
		
		for(String key : map.keySet()){
			String filename = args[1]+"_"+(key.replaceAll(",", "_"))+".fasta";
			DataSet ds = new DataSet("",map.get(key));
			ds.save(new FileOutputStream(filename), '>', new SimpleSequenceAnnotationParser());
		}
		
	}

}
