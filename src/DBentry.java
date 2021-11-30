/**
 * Type of entry of the databse
 */
public class DBentry {
    //jai mi la sequence dans 1 une seule string, mais on verra apres si sa serait mieux des la mettre dans une
    //arraylist caractere par caracteres
    String sequenceInfo;
    String sequence;

    public DBentry(String sequenceInfo, String sequence){
        this.sequenceInfo = sequenceInfo;
        this.sequence = sequence;
    }

    public String getSequence() {
        return sequence;
    }

    public String getSequenceInfo() {
        return sequenceInfo;
    }
}
