public class HspStat {


    String currentHsp;
    double lambda = 0.192;
    double k = 0.176;
    int bruteScore;
    double bitScore;
    double eValue;
    String database;
    int positionDatabase;
    int m;
    int n;
    int positionInput;

    public HspStat(String input,String currentHsp, int positionInput, int positionDatabase, String database, int m){
        this.m = m;
        this.n = input.length();
        this.currentHsp = currentHsp;
        this.database = database;
        this.positionInput = positionInput;
        this.positionDatabase = positionDatabase;

        this.bruteScore();
        this.bitScore = Math.round((this.lambda * this.bruteScore - Math.log(k))/Math.log(2));
        this.eValue = this.m * this.n / Math.pow(2, this.bitScore);
    }

    public void bruteScore(){
        this.bruteScore = 0;
        for (int i = 0; i < this.currentHsp.length(); i++){
            if (this.currentHsp.charAt(i) == this.database.charAt(i+positionDatabase)){
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
}
