package projects.sigma;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedList;

import de.jstacs.data.DNADataSet;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;

public class Filter {

	public static void main(String[] args) throws FileNotFoundException, WrongAlphabetException, EmptyDataSetException, WrongLengthException, IOException {
		DNADataSet data = new DNADataSet(args[0], '>', new SimpleSequenceAnnotationParser());
		
		HashSet<String> set = new HashSet<>();
		
		BufferedReader read = new BufferedReader(new FileReader(args[1]));
		
		String str = null;
		while( (str = read.readLine()) != null ){
			set.add(str);
		}
		read.close();
		
		LinkedList<Sequence> in = new LinkedList<>();
		LinkedList<Sequence> out = new LinkedList<>();
		
		for(int i=0;i<data.getNumberOfElements();i++){
			Sequence seq = data.getElementAt(i);
			String ann = (String) seq.getAnnotation()[0].getResultForName("unparsed comment").getValue();
			if(set.contains(ann)){
				in.add(seq);
			}else{
				out.add(seq);
			}
		}
		
		
		DataSet ds = new DataSet("",in);
		ds.save(new FileOutputStream(args[0]+"_in.fa"), '>', new SimpleSequenceAnnotationParser());
		
		ds = new DataSet("",out);
		ds.save(new FileOutputStream(args[0]+"_out.fa"), '>', new SimpleSequenceAnnotationParser());
		

	}

}
