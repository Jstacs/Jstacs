package projects.xanthogenomes;

import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;

public class FindInvertedRepeat {

	public static void main(String[] args) throws Exception{
		
		Sequence seq = Sequence.create(DNAAlphabetContainer.SINGLETON, args[0]);
		
		for(int i=0;i<seq.getLength();i++){
			int minK = 4;
			for(int j=i+1;j<seq.getLength();j++){
				if(i==0||j==seq.getLength()-1||seq.discreteVal(i-1) != 3-seq.discreteVal(j+1)){
					int k=0;			
					while( i+k<seq.getLength() && seq.discreteVal(i+k) == 3-seq.discreteVal(j-k) && i+k<j-k ){
						k++;
					}
					if(k > minK){
						System.out.println(i+" "+j+" "+(j-i+1)+" "+k+" "+seq.getSubSequence(i, k));
					}
				}
			}
			
		}
		
		
	}

}
