package com.github.micycle1.wevo4j;

import java.util.*;
import org.locationtech.jts.geom.*;

import java.util.*;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.operation.polygonize.Polygonizer;
import org.locationtech.jts.operation.union.UnaryUnionOp;

public final class MwvdExporter {
	private final MwvdGeom g;
	private final GeometryFactory gf;
	private final double chordTol;

	public MwvdExporter(MwvdGeom g, GeometryFactory gf, double chordTol) {
		this.g = g;
		this.gf = gf;
		this.chordTol = chordTol;
	}

	public sealed interface EdgePiece permits ArcPiece, SegPiece {
	}

	public record ArcPiece(MwvdGeom.Arc arc, int s1, int s2) implements EdgePiece {
	}

	public record SegPiece(MwvdGeom.Segment seg, int s1, int s2) implements EdgePiece {
	}

	public GeometryCollection exportRegions(List<EdgePiece> edges, List<MwvdModel.Site> sites, Envelope clip) {
		Polygon clipPoly = (Polygon) gf.toGeometry(clip);
		List<Geometry> linework = new ArrayList<>(edges.size() + 1);

		for (EdgePiece e : edges) {
			Geometry curve = toGeometry(e);
			Geometry clipped = curve.intersection(clipPoly);
			collectLinear(clipped, linework);
		}

		collectLinear(clipPoly.getBoundary(), linework);

		Geometry noded = UnaryUnionOp.union(linework);

		Polygonizer polygonizer = new Polygonizer();
		polygonizer.add(noded);

		@SuppressWarnings("unchecked")
		Collection<Polygon> faces = polygonizer.getPolygons();

		// disabled temporarily
//		Map<Integer, List<Geometry>> perSite = new TreeMap<>();
//		for (MwvdModel.Site s : sites) {
//			perSite.put(s.id, new ArrayList<>());
//		}
//
//		for (Polygon face : faces) {
//			if (face.isEmpty())
//				continue;
//			Coordinate q = face.getInteriorPoint().getCoordinate();
//			MwvdModel.Site owner = owner(q, sites);
//			perSite.get(owner.id).add(face);
//		}
//
//		List<Geometry> out = new ArrayList<>(sites.size());
//		for (MwvdModel.Site s : sites) {
//			Geometry cell = perSite.get(s.id).isEmpty() ? gf.createPolygon() : UnaryUnionOp.union(perSite.get(s.id));
//			cell = cell.intersection(clipPoly);
//			cell.setUserData(s.id);
//			out.add(cell);
//		}

		return gf.createGeometryCollection(linework.toArray(Geometry[]::new));
	}

	private MwvdModel.Site owner(Coordinate q, List<MwvdModel.Site> sites) {
		MwvdModel.Site best = null;
		double bestT = Double.POSITIVE_INFINITY;
		for (MwvdModel.Site s : sites) {
			double t = s.sqTimeAt(q);
			if (best == null || g.lt(t, bestT) || (g.eq(t, bestT) && s.id < best.id)) {
				best = s;
				bestT = t;
			}
		}
		return best;
	}

	private Geometry toGeometry(EdgePiece e) {
		return switch (e) {
			case ArcPiece a -> gf.createLineString(approxArc(a.arc()));
			case SegPiece s -> gf.createLineString(new Coordinate[] { s.seg().a(), s.seg().b() });
		};
	}

	private void collectLinear(Geometry g0, List<Geometry> out) {
		if (g0 == null || g0.isEmpty())
			return;

		if (g0 instanceof LineString || g0 instanceof MultiLineString) {
			out.add(g0);
			return;
		}

		for (int i = 0; i < g0.getNumGeometries(); i++) {
			collectLinear(g0.getGeometryN(i), out);
		}
	}

	private Coordinate[] approxArc(MwvdGeom.Arc arc) {
		double sweep = g.ccwSweep(arc.a0(), arc.a1());
		if (g.eq(sweep, 0.0)) {
			return new Coordinate[] { g.pointOnCircle(arc.circle(), arc.a0()), g.pointOnCircle(arc.circle(), arc.a1()) };
		}

		double r = arc.circle().r();
		double maxStep = r <= chordTol ? Math.PI / 8.0 : Math.min(Math.PI / 8.0, 2.0 * Math.acos(Math.max(-1.0, 1.0 - chordTol / r)));

		int n = Math.max(1, (int) Math.ceil(sweep / maxStep));
		Coordinate[] cs = new Coordinate[n + 1];

		for (int i = 0; i <= n; i++) {
			double a = g.normAngle(arc.a0() + sweep * i / n);
			cs[i] = g.pointOnCircle(arc.circle(), a);
		}
		return cs;
	}
}