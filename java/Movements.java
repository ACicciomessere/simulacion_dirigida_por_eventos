public class Movements {
    private Particle[] particles;

    public Movements(Particle[] particles){
        this.particles = particles;
    }

    public void update(double dt, double L, double r0, boolean periodic){
        for (Particle p : particles) {
            p.updatePosition(L, r0, periodic);
        }
    }

    public void updateTheta(double theta){
        for (Particle p : particles) {
            p.setTheta(theta);
        }
    }

    
}