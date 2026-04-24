import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.PriorityQueue;

public class Collisions {
    private PriorityQueue<Event> pq;
    private double simTime;
    private final List<Particle> particles;

    private final double R_OUTER;
    private final double R_INNER;
    private final double CX;
    private final double CY;

    public Collisions(List<Particle> particles, double rOuter, double rInner, double cx, double cy) {
        this.particles = particles;
        this.R_OUTER = rOuter;
        this.R_INNER = rInner;
        this.CX = cx;
        this.CY = cy;
        this.simTime = 0.0;
        this.pq = new PriorityQueue<>();
    }

    private void predict(Particle p) {
        if (p == null) return;

        double tOuter = p.timeToHitCircle(CX, CY, R_OUTER, false);
        if (tOuter > 1e-12)
            pq.add(new Event(simTime + tOuter, p, null));

        double tInner = p.timeToHitCircle(CX, CY, R_INNER, true);
        if (tInner > 1e-12)
            pq.add(new Event(simTime + tInner, p, p));
    }

    public void simulate(double maxTime, int nEventsPerPrint, String outputPath) throws IOException {

        for (Particle p : particles)
            predict(p);

        pq.add(new Event(0, null, null));

        File outFile = new File(outputPath);
        File parent = outFile.getParentFile();
        if (parent != null && !parent.exists())
            parent.mkdirs();

        int N = particles.size();
        // Pre-size: time header (12) + N particles (~30 chars each)
        StringBuilder sb = new StringBuilder(12 + N * 30);

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outFile), 1 << 20)) {
            bw.write("# t\n# x y vx vy\n");

            int countEvents = 0;

            while (!pq.isEmpty()) {
                Event event = pq.poll();
                if (!event.isValid()) continue;
                if (event.time < simTime - 1e-12) continue;

                // Lazy update: only move the 1-2 particles involved in this event
                Particle a = event.a;
                Particle b = event.b;
                if (a != null) a.moveTo(event.time);
                if (b != null && b != a) b.moveTo(event.time);

                simTime = event.time;
                countEvents++;

                if (countEvents % nEventsPerPrint == 0) {
                    sb.setLength(0);
                    appendDouble6f(sb, simTime); sb.append('\n');
                    for (Particle p : particles) {
                        appendDouble6f(sb, p.getX(simTime)); sb.append(' ');
                        appendDouble6f(sb, p.getY(simTime)); sb.append(' ');
                        appendDouble6f(sb, p.vx);            sb.append(' ');
                        appendDouble6f(sb, p.vy);            sb.append('\n');
                    }
                    bw.write(sb.toString());
                }

                if (simTime > maxTime) break;

                // Resolve collision
                if (a != null && b != null) {
                    if (a == b)
                        a.bounceOffCircle(CX, CY);
                    else
                        a.bounceOff(b);
                } else if (a != null) {
                    a.bounceOffCircle(CX, CY);
                }

                predict(a);
                if (b != null && b != a) predict(b);
            }
        }

        System.out.println("Simulation finished.");
    }

    // Fast "%.6f" formatter using integer arithmetic — ~10x faster than String.format
    private static void appendDouble6f(StringBuilder sb, double val) {
        if (val < 0) { sb.append('-'); val = -val; }
        long scaled = Math.round(val * 1_000_000.0);
        long intPart = scaled / 1_000_000L;
        int  fracPart = (int)(scaled % 1_000_000L);
        sb.append(intPart).append('.');
        if (fracPart < 100000) sb.append('0');
        if (fracPart < 10000)  sb.append('0');
        if (fracPart < 1000)   sb.append('0');
        if (fracPart < 100)    sb.append('0');
        if (fracPart < 10)     sb.append('0');
        sb.append(fracPart);
    }
}
