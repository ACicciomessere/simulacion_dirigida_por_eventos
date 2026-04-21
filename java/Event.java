
public class Event implements Comparable<Event> {
    public final double time;
    public final Particle a, b;
    private final int countA, countB;

    public Event(double time, Particle a, Particle b) {
        this.time = time;
        this.a = a;
        this.b = b;
        this.countA = (a != null) ? a.collisionCount : -1;
        this.countB = (b != null) ? b.collisionCount : -1;
    }

    @Override
    public int compareTo(Event that) {
        return Double.compare(this.time, that.time);
    }

    public boolean isValid() {
        if (a != null && a.collisionCount != countA) return false;
        return b == null || b.collisionCount == countB;
    }
}



