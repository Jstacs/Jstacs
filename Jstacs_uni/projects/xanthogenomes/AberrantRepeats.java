package projects.xanthogenomes;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import projects.xanthogenomes.TALE.Repeat;
import projects.xanthogenomes.TALE.Type;

public class AberrantRepeats {

	public static void main(String[] args) throws Exception {
		
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder(FileManager.readFile(args[0]));

		TALE[] tales = builder.getAllTALEs();
		
		for(int i=0;i<tales.length;i++){
			TALE dna = tales[i].getDnaOriginal();
			for(int j=0;j<dna.getNumberOfRepeats()-1;j++){
				Repeat r = dna.getRepeat(j);
				
				int pos = tales[i].getRepeat(j).getRvdPosition();
				int len = tales[i].getRepeat(j).getRvdLength();
				
				if(tales[i].getRepeat(j).getType()==Type.SHORT){
					System.out.println((">"+dna.getId()+" repeat "+j).replaceAll(" ", "_"));
					Sequence seq = r.getRepeat();
					System.out.println(seq.getSubSequence(0, pos*3)+"NNN"+seq.getSubSequence((pos+len)*3));
					//System.out.println(seq);
					//System.out.println(pos+" "+len);
				}
			}
		}
		
	}

}
