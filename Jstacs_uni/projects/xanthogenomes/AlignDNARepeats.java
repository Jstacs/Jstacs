package projects.xanthogenomes;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.FileManager;
import de.jstacs.utils.IntList;

public class AlignDNARepeats {

	public static void main(String[] args) throws Exception {
		
		//TalCS1 vs TalBB
		//String first = "TalCS1";
		//String[] second = new String[]{"TalBB3","TalBB1","TalBB2", "TalBB7"};
		
		//TalDS1 vs TalAC
		//String first = "TalDS1";
		//String[] second = new String[]{"TalAC9","TalAC6","TalAC8","TalAC3","TalAC1","TalAC7","TalAC4","TalAC5"};
		
		//TalAI vs TalDR
		//String first = "TalAI2";
		//String[] second = new String[]{"TalDR1","TalDR2"};
		
		//TalDU1 vs TalAQ
		//String first = "TalDU1";
		//String[] second = new String[]{"TalAQ5","TalAQ11","TalAQ9"};
		
		//TalDT1 vs TalAR
		//String first = "TalDT1";
		//String[] second = new String[]{"TalAR3","TalAR4","TalAR6"};
		
		//TalCM1 vs TalCM
		String first = "TalCM5";
		String[] second = new String[]{"TalCM2","TalCM4"};
		
		TALEFamilyBuilder builder = new TALEFamilyBuilder(FileManager.readFile(args[0]));
		
		TALE ft = null;
		TALE[] st = new TALE[second.length];
		TALE[] tales = builder.getAllTALEs();
		for(int i=0;i<tales.length;i++){
			if(tales[i].getId().startsWith(first)){
				ft = tales[i];
			}else{
				for(int j=0;j<second.length;j++){
					if(tales[i].getId().startsWith(second[j])){
						st[j] = tales[i];
					}
				}
			}
		}
		
		IntList firstUsed = new IntList();
		//TalCS1 vs TalBB
		/*for(int i=0;i<4;i++){
			firstUsed.add(i);
		}
		for(int i=14;i<ft.getNumberOfRepeats();i++){
			firstUsed.add(i);
		}*/
		
		//TalAI vs TalDR
		/*for(int i=0;i<12;i++){
			firstUsed.add(i);
		}
		for(int i=13;i<ft.getNumberOfRepeats();i++){
			firstUsed.add(i);
		}*/
		
		//TalCM1 vs TalCM
		for(int i=0;i<3;i++){
			firstUsed.add(i);
		}
		
		for(int i=4;i<11;i++){
			firstUsed.add(i);
		}
		
		for(int i=12;i<ft.getNumberOfRepeats();i++){
			firstUsed.add(i);
		}
		
		/*for(int i=0;i<ft.getNumberOfRepeats();i++){
			firstUsed.add(i);
		}*/
		
		
		IntList secondUsed = new IntList();
		//TalDS1 vs TalAC
		/*for(int i=0;i<6;i++){
			secondUsed.add(i);
		}
		for(int i=7;i<st[0].getNumberOfRepeats();i++){
			secondUsed.add(i);
		}*/
		
		//TalDU1 vs TalAQ
		/*for(int i=0;i<12;i++){
			secondUsed.add(i);
		}
		for(int i=13;i<st[0].getNumberOfRepeats();i++){
			secondUsed.add(i);
		}*/
		
		//TalDT1 vs TalAR
		/*for(int i=0;i<8;i++){
			secondUsed.add(i);
		}
		for(int i=9;i<st[0].getNumberOfRepeats();i++){
			secondUsed.add(i);
		}*/
		for(int i=0;i<st[0].getNumberOfRepeats();i++){
			secondUsed.add(i);
		}
		
		if(firstUsed.length() != secondUsed.length()){
			throw new RuntimeException(firstUsed.length()+" "+secondUsed.length());
		}
		
		Sequence fn = ft.getDnaOriginal().getStart();
		Sequence[] sns = new Sequence[st.length];
		for(int i=0;i<st.length;i++){
			sns[i] = st[i].getDnaOriginal().getStart();
		}
		
		int offset = 0;
		compare(fn,sns,offset,"n",false);
		offset += fn.getLength();
		
		
		for(int i=0;i<firstUsed.length();i++){
			Sequence fr = ft.getDnaOriginal().getRepeat(firstUsed.get(i)).getRepeat();
			Sequence[] srs = new Sequence[st.length];
			for(int j=0;j<st.length;j++){
				srs[j] = st[j].getDnaOriginal().getRepeat(secondUsed.get(i)).getRepeat();
			}
			compare(fr,srs,offset,"r"+i,i>0 && (firstUsed.get(i)-1!=firstUsed.get(i-1) || secondUsed.get(i)-1 != secondUsed.get(i-1)));
			offset += fr.getLength();
		}
		
		
		
		Sequence fc = ft.getDnaOriginal().getEnd();
		Sequence[] scs = new Sequence[st.length];
		for(int i=0;i<st.length;i++){
			scs[i] = st[i].getDnaOriginal().getEnd();
		}
		compare(fc,scs,offset,"c",false);
		
	}

	private static void compare(Sequence fn, Sequence[] sns, int offset,String type, boolean breakpoint) {
		int[] diffs = new int[fn.getLength()];
		for(int i=0;i<fn.getLength();i++){
			for(int j=0;j<sns.length;j++){
				if(fn.discreteVal(i) != sns[j].discreteVal(i)){
					diffs[i]++;
				}
			}
		}
		for(int i=0;i<diffs.length;i++){
			System.out.println((i+offset)+"\t"+diffs[i]+"\t"+type+"\t"+(i==0 ? "y" : "n")+"\t"+(i==0 && breakpoint ? "b":"n"));
		}
	}

}
