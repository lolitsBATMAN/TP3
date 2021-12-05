public class HspStat {


    private final String currentHsp;
    private int bruteScore;
    private final double bitScore;
    private final double eValue;
    private final DBentry database;
    private final int positionDatabase;
    private int m;
    private int n;
    private int positionInput;

    public HspStat(String input,String currentHsp, int positionInput, int positionDatabase, DBentry database, int m){
        this.m = m;
        this.n = input.length();
        this.currentHsp = currentHsp;
        this.database = database;
        this.positionInput = positionInput;
        this.positionDatabase = positionDatabase;

        this.bruteScore();
        this.bitScore = Math.round((0.192 * this.bruteScore - Math.log(0.176))/Math.log(2));
        this.eValue = this.m * this.n / Math.pow(2, this.bitScore);
    }

    public void bruteScore(){
        this.bruteScore = 0;
        for (int i = 0; i < this.currentHsp.length(); i++){
            if (this.currentHsp.charAt(i) == this.database.getSequence().charAt(i+positionDatabase)){
                this.bruteScore += 5;
            } else {
                this.bruteScore -= 4;
            }
        }
    }

    public double getBruteScore(){
        return this.bruteScore;
    }

    public double geteValue(){
        return this.eValue;
    }

    public double getBitScore(){
        return this.bitScore;
    }

    public String getCurrentHsp(){return this.currentHsp;}

    public DBentry getDatabase(){return this.database;}

    public int getPositionDatabase(){return this.positionDatabase;}

    public int getPositionInput(){return this.positionInput;}
}
