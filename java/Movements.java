public class Movements {
    private Particle[] particles;

    public Movements(Particle[] particles){
        this.particles = particles;
    }

    public void update(double dt, double r0, boolean periodic){
        for (Particle p : particles) {
            p.updatePosition(r0, 1.0, periodic);
        }
    }

    public void updateTheta(double theta){
        for (Particle p : particles) {
            p.setTheta(theta);
        }
    }

    
}