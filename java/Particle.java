public class Particle {
    public double x, y, vx, vy, radius, mass;
    public int collisionCount;

    public Particle(double x, double y, double vx, double vy, double radius, double mass) {
        this.x = x;
        this.y = y;
        this.vx = vx;
        this.vy = vy;
        this.radius = radius;
        this.mass = mass;
        this.collisionCount = 0;
    }

    public void move(double dt) {
        this.x += this.vx * dt;
        this.y += this.vy * dt;
    }

    public double timeToHit(Particle other) {
        if (this == other)
            return Double.POSITIVE_INFINITY;

        double dx = other.x - this.x;
        double dy = other.y - this.y;
        double dvx = other.vx - this.vx;
        double dvy = other.vy - this.vy;

        double dvdr = dvx * dx + dvy * dy;
        if (dvdr > 0)
            return Double.POSITIVE_INFINITY;

        double dvdv = dvx * dvx + dvy * dvy;
        if (dvdv == 0)
            return Double.POSITIVE_INFINITY;

        double drdr = dx * dx + dy * dy;
        double sigma = this.radius + other.radius;

        double d = (dvdr * dvdr) - dvdv * (drdr - sigma * sigma);
        if (d < 0)
            return Double.POSITIVE_INFINITY;

        double t = -(dvdr + Math.sqrt(d)) / dvdv;
        return t > 0 ? t : Double.POSITIVE_INFINITY;
    }

    public double timeToHitVerticalWall(double xMin, double xMax) {
        if (vx > 0)
            return (xMax - radius - x) / vx;
        if (vx < 0)
            return (xMin + radius - x) / vx;
        return Double.POSITIVE_INFINITY;
    }

    public double timeToHitHorizontalWall(double yMin, double yMax) {
        if (vy > 0)
            return (yMax - radius - y) / vy;
        if (vy < 0)
            return (yMin + radius - y) / vy;
        return Double.POSITIVE_INFINITY;
    }

    public double timeToHitPoint(double px, double py) {
        double dx = px - this.x;
        double dy = py - this.y;
        double dvx = this.vx;
        double dvy = this.vy;

        double dvdr = dx * dvx + dy * dvy;
        if (dvdr <= 0)
            return Double.POSITIVE_INFINITY;

        double dvdv = dvx * dvx + dvy * dvy;
        if (dvdv == 0)
            return Double.POSITIVE_INFINITY;

        double drdr = dx * dx + dy * dy;
        double sigma = this.radius;

        double d = (dvdr * dvdr) - dvdv * (drdr - sigma * sigma);
        if (d < 0)
            return Double.POSITIVE_INFINITY;

        double t = -(dvdr + Math.sqrt(d)) / dvdv;
        return t > 0 ? t : Double.POSITIVE_INFINITY;
    }

    public void bounceOff(Particle other) {
        double dx = other.x - this.x;
        double dy = other.y - this.y;
        double dvx = other.vx - this.vx;
        double dvy = other.vy - this.vy;

        double dvdr = dx * dvx + dy * dvy;
        double dist = this.radius + other.radius;

        double impulse = (2 * this.mass * other.mass * dvdr) / ((this.mass + other.mass) * dist);
        double Jx = impulse * dx / dist;
        double Jy = impulse * dy / dist;

        this.vx += Jx / this.mass;
        this.vy += Jy / this.mass;
        other.vx -= Jx / other.mass;
        other.vy -= Jy / other.mass;

        this.collisionCount++;
        other.collisionCount++;
    }

    public void bounceOffCircle(double cx, double cy) {
        double dx = this.x - cx;
        double dy = this.y - cy;
        double dist = Math.sqrt(dx * dx + dy * dy);

        double nx = dx / dist;
        double ny = dy / dist;

        double vdotn = vx * nx + vy * ny;

        vx -= 2 * vdotn * nx;
        vy -= 2 * vdotn * ny;

        collisionCount++;
    }

    // cuánto tiempo falta para la colisión
    public double timeToHitCircle(double cx, double cy, double R, boolean isInnerWall) {
        double dx = this.x - cx;
        double dy = this.y - cy;
        double dvx = this.vx;
        double dvy = this.vy;

        double targetRadius = isInnerWall ? R + this.radius : R - this.radius;
        double a = dvx * dvx + dvy * dvy;
        double b = 2 * (dx * dvx + dy * dvy);
        double c = dx * dx + dy * dy - targetRadius * targetRadius;

        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0)
            return Double.POSITIVE_INFINITY;

        double t1 = (-b - Math.sqrt(discriminant)) / (2 * a);
        double t2 = (-b + Math.sqrt(discriminant)) / (2 * a);

        if (!isInnerWall) {
            // Partícula DENTRO del recinto: siempre usar la raíz MAYOR (t2)
            // t1 es la intersección "de entrada" que ya quedó atrás
            return t2 > 0 ? t2 : Double.POSITIVE_INFINITY;
        } else {
            // Partícula FUERA del obstáculo: usar la raíz MENOR (t1)
            if (t1 > 0)
                return t1;
            if (t2 > 0)
                return t2;
            return Double.POSITIVE_INFINITY;
        }
    }

}
