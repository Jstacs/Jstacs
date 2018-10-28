package projects.xanthogenomes;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.parameters.FileParameter;
import projects.xanthogenomes.TALEFamilyBuilder.TALEFamily;

public class Frameshifts {

	public static void main(String[] args) throws Exception {
		
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder(FileManager.readFile(args[0]));

		TALE[] tales = builder.getAllTALEs();
		
		for(int i=0;i<tales.length-1;i++){
			for(int j=i+1;j<tales.length;j++){
				Sequence rvds1 = tales[i].getRvdSequence();
				Sequence rvds2 = tales[j].getRvdSequence();
				if(rvds1.getLength() != rvds2.getLength()
						&& rvds1.getSubSequence(0, 4).getHammingDistance(rvds2.getSubSequence(0, 4)) == 0
						&& rvds1.getSubSequence(rvds1.getLength()-4).getHammingDistance(rvds2.getSubSequence(rvds2.getLength()-4))==0){
					System.out.println(tales[i].getId());
					System.out.println(rvds1);
					System.out.println(tales[j].getId());
					System.out.println(rvds2);
					System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++");
				}
				
			}
		}
		
		
		System.out.println("#######################################################");
		
		
		for(int i=0;i<tales.length-1;i++){
			for(int j=i+1;j<tales.length;j++){
				Sequence rvds1 = tales[i].getRvdSequence();
				Sequence rvds2 = tales[j].getRvdSequence();
				
				if(rvds1.getLength()<rvds2.getLength()){
					Sequence temp = rvds1;
					rvds1 = rvds2;
					rvds2 = temp;
				}
				
				if( rvds1.getLength() > rvds2.getLength() ){
					if(rvds1.getSubSequence(0,rvds2.getLength()).getHammingDistance(rvds2) == 0 ||
							rvds1.getSubSequence(rvds1.getLength()-rvds2.getLength(), rvds2.getLength()).getHammingDistance(rvds2) == 0){
						System.out.println(tales[i].getId());
						System.out.println(rvds1);
						System.out.println(tales[j].getId());
						System.out.println(rvds2);
						System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++");
					}
				}
				
				
			}
		}
		
		
		TALEFamily[] fams = builder.getFamilies();
		for(int i=0;i<fams.length;i++){
			String fam = fams[i].getFamilyId();
			TALE[] ts = fams[i].getFamilyMembers();
			for(int j=0;j<ts.length;j++){
				if(ts[j].getStrain().startsWith("Xoc")){
					System.out.println(fam+"\t"+ts[j].getStrain()+"\t"+ts[j].getId());
				}
			}
		}
		
		
		
	}

}
