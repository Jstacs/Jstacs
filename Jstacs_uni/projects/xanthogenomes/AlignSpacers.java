package projects.xanthogenomes;

import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.cost.SimpleCosts;
import de.jstacs.data.DNADataSet;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.annotation.SimpleSequenceAnnotationParser;

public class AlignSpacers {

	public static void main(String[] args) throws Exception {
		
		DNADataSet ds = new DNADataSet(args[0], '>', new SimpleSequenceAnnotationParser());
		
		System.out.print("name");
		for(int i=0;i<ds.getNumberOfElements();i++){
			System.out.print("\t"+ds.getElementAt(i).getAnnotation()[0].getResultAt(0).getValue());
		}
		System.out.println();
		
		Alignment al = new Alignment(new SimpleCosts(0, 1, 1));
		
		for(int i=0;i<ds.getNumberOfElements();i++){
			System.out.print(ds.getElementAt(i).getAnnotation()[0].getResultAt(0).getValue());
			for(int j=0;j<ds.getNumberOfElements();j++){
				Sequence s1 = ds.getElementAt(i);
				Sequence s2 = ds.getElementAt(j);
				System.out.print("\t"+(al.getAlignment(AlignmentType.GLOBAL, s1, s2).getCost()));
			}
			System.out.println();
		}

	}

}
