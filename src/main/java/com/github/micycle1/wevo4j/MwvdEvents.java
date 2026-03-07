package com.github.micycle1.wevo4j;

import java.util.Comparator;
import java.util.PriorityQueue;
import org.locationtech.jts.geom.Coordinate;

public final class MwvdEvents {
    private MwvdEvents() {}

    public enum Type { COLL, EDGE, DOM }

    public abstract static class Event {
        public final Coordinate p;
        public final double sqTime;
        public final MwvdModel.Site site;

        protected Event(Coordinate p, double sqTime, MwvdModel.Site site) {
            this.p = new Coordinate(p.x, p.y);
            this.sqTime = sqTime;
            this.site = site;
        }

        public abstract Type type();
    }

    public static final class CollEvent extends Event {
        public final MwvdModel.Trajectory traj1;
        public final MwvdModel.Trajectory traj2;
        public final boolean pierces;

        public CollEvent(Coordinate p, double t, MwvdModel.Site site,
                         MwvdModel.Trajectory traj1, MwvdModel.Trajectory traj2, boolean pierces) {
            super(p, t, site);
            this.traj1 = traj1; this.traj2 = traj2; this.pierces = pierces;
        }
        @Override public Type type() { return Type.COLL; }
    }

    public static final class EdgeEvent extends Event {
        public final MwvdModel.MovingIntersection isect1;
        public final MwvdModel.MovingIntersection isect2;

        public EdgeEvent(Coordinate p, double t, MwvdModel.Site site,
                         MwvdModel.MovingIntersection isect1, MwvdModel.MovingIntersection isect2) {
            super(p, t, site);
            this.isect1 = isect1; this.isect2 = isect2;
        }
        @Override public Type type() { return Type.EDGE; }
    }

    public static final class DomEvent extends Event {
        public final MwvdModel.MovingIntersection isect1;
        public final MwvdModel.MovingIntersection isect2;

        public DomEvent(Coordinate p, double t, MwvdModel.Site site,
                        MwvdModel.MovingIntersection isect1, MwvdModel.MovingIntersection isect2) {
            super(p, t, site);
            this.isect1 = isect1; this.isect2 = isect2;
        }
        @Override public Type type() { return Type.DOM; }
    }

    public static PriorityQueue<Event> newQueue(MwvdGeom g) {
        return new PriorityQueue<>(new EventComparator(g));
    }

    /** Deterministic total order using EPS buckets then stable raw/id ties. */
    public static final class EventComparator implements Comparator<Event> {
        private final MwvdGeom g;
        public EventComparator(MwvdGeom g) { this.g = g; }

        @Override public int compare(Event a, Event b) {
            int c = compareScalar(a.sqTime, b.sqTime);
            if (c != 0) return c;
            c = compareScalar(a.p.x, b.p.x);
            if (c != 0) return c;
            c = compareScalar(a.p.y, b.p.y);
            if (c != 0) return c;
            c = Integer.compare(a.type().ordinal(), b.type().ordinal());
            if (c != 0) return c;
            c = Integer.compare(a.site.id, b.site.id);
            if (c != 0) return c;
            return payloadCompare(a, b);
        }

        private int compareScalar(double x, double y) {
            if (!g.eq(x, y)) return Double.compare(x, y);
            return 0;
        }

        private int payloadCompare(Event a, Event b) {
            if (a instanceof CollEvent ca && b instanceof CollEvent cb) {
                int c = idCmp(ca.traj1.id(), cb.traj1.id());
                if (c != 0) return c;
                c = Boolean.compare(ca.pierces, cb.pierces);
                if (c != 0) return c;
                return Boolean.compare(ca.traj1.isLeft, cb.traj1.isLeft);
            }
            if (a instanceof DomEvent da && b instanceof DomEvent db) {
                return isectCmp(da.isect1.id(), db.isect1.id());
            }
            if (a instanceof EdgeEvent ea && b instanceof EdgeEvent eb) {
                int c = isectCmp(ea.isect1.id(), eb.isect1.id());
                if (c != 0) return c;
                return isectCmp(ea.isect2.id(), eb.isect2.id());
            }
            return 0;
        }

        private int idCmp(MwvdModel.SitePair a, MwvdModel.SitePair b) {
            int c = Integer.compare(a.a(), b.a());
            return c != 0 ? c : Integer.compare(a.b(), b.b());
        }

        private int isectCmp(MwvdModel.MovIsectId a, MwvdModel.MovIsectId b) {
            int c = Integer.compare(a.s1(), b.s1());
            if (c != 0) return c;
            c = Integer.compare(a.s2(), b.s2());
            if (c != 0) return c;
            c = Boolean.compare(a.isLeft(), b.isLeft());
            if (c != 0) return c;
            return Boolean.compare(a.isFirst(), b.isFirst());
        }
    }
}
