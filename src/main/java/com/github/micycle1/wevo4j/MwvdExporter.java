package com.github.micycle1.wevo4j;

import java.util.*;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Polygon;
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
		Coordinate start();
		Coordinate end();
		int s1();
		int s2();
	}

	public record ArcPiece(MwvdGeom.Arc arc, Coordinate start, Coordinate end, int s1, int s2) implements EdgePiece {
		public ArcPiece {
			start = new Coordinate(start.x, start.y);
			end = new Coordinate(end.x, end.y);
		}
	}

	public record SegPiece(MwvdGeom.Segment seg, Coordinate start, Coordinate end, int s1, int s2) implements EdgePiece {
		public SegPiece {
			start = new Coordinate(start.x, start.y);
			end = new Coordinate(end.x, end.y);
		}
	}

	public GeometryCollection exportRegions(List<EdgePiece> edges, List<MwvdModel.Site> sites, Envelope clip) {
		Polygon clipPoly = (Polygon) gf.toGeometry(clip);
		List<EdgePiece> clipped = new ArrayList<>(edges.size());
		List<Coordinate> boundaryVertices = new ArrayList<>();

		for (EdgePiece edge : edges) {
			for (EdgePiece piece : clipToBox(edge, clip)) {
				clipped.add(piece);
				if (onBoundary(piece.start(), clip)) {
					addUnique(boundaryVertices, piece.start());
				}
				if (onBoundary(piece.end(), clip)) {
					addUnique(boundaryVertices, piece.end());
				}
			}
		}

		for (Coordinate corner : boxCorners(clip)) {
			addUnique(boundaryVertices, corner);
		}

		HalfEdgeGraph graph = new HalfEdgeGraph();
		for (EdgePiece piece : clipped) {
			graph.add(piece);
		}
		addBoundaryEdges(graph, boundaryVertices, clip);

		List<Polygon> faces = graph.walkFaces();
		Map<Integer, List<Geometry>> perSite = new TreeMap<>();
		for (MwvdModel.Site s : sites) {
			perSite.put(s.id, new ArrayList<>());
		}

		for (Polygon face : faces) {
			if (face.isEmpty()) {
				continue;
			}
			Coordinate q = face.getInteriorPoint().getCoordinate();
			MwvdModel.Site owner = owner(q, sites);
			perSite.get(owner.id).add(face);
		}

		List<Geometry> out = new ArrayList<>(sites.size());
		for (MwvdModel.Site s : sites) {
			Geometry cell = perSite.get(s.id).isEmpty() ? gf.createPolygon() : UnaryUnionOp.union(perSite.get(s.id));
			cell = cell.intersection(clipPoly);
			out.add(cell);
		}

		carveContainedCells(out);
		for (int i = 0; i < sites.size(); i++) {
			out.get(i).setUserData(sites.get(i).id);
		}

		return gf.createGeometryCollection(out.toArray(Geometry[]::new));
	}

	private void carveContainedCells(List<Geometry> cells) {
		List<Integer> order = new ArrayList<>(cells.size());
		for (int i = 0; i < cells.size(); i++) {
			order.add(i);
		}
		order.sort(Comparator.comparingDouble((Integer i) -> cells.get(i).getArea()));

		for (int childIdx : order) {
			Geometry child = cells.get(childIdx);
			if (child.isEmpty()) {
				continue;
			}

			Geometry probe = child.getInteriorPoint();
			for (int parentIdx : order) {
				if (parentIdx == childIdx) {
					continue;
				}

				Geometry parent = cells.get(parentIdx);
				if (parent.isEmpty() || parent.getArea() <= child.getArea()) {
					continue;
				}
				if (!parent.covers(probe)) {
					continue;
				}

				cells.set(parentIdx, parent.difference(child));
			}
		}
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

	private List<EdgePiece> clipToBox(EdgePiece edge, Envelope box) {
		return switch (edge) {
			case ArcPiece a -> clipArc(a, box);
			case SegPiece s -> clipSegment(s, box);
		};
	}

	private List<EdgePiece> clipSegment(SegPiece edge, Envelope box) {
		List<Coordinate> cuts = new ArrayList<>();
		addUnique(cuts, edge.start());
		addUnique(cuts, edge.end());

		for (MwvdGeom.Segment side : boxSides(box)) {
			for (Coordinate p : g.segmentSegmentIntersections(edge.seg(), side)) {
				if (g.pointOnSegment(p, edge.seg())) {
					addUnique(cuts, p);
				}
			}
		}

		cuts.sort(Comparator.comparingDouble(p -> segmentParam(edge.seg(), p)));
		List<EdgePiece> out = new ArrayList<>();
		for (int i = 0; i + 1 < cuts.size(); i++) {
			Coordinate a = cuts.get(i);
			Coordinate b = cuts.get(i + 1);
			if (g.samePoint(a, b)) {
				continue;
			}
			double ta = segmentParam(edge.seg(), a);
			double tb = segmentParam(edge.seg(), b);
			Coordinate mid = segmentPoint(edge.seg(), 0.5 * (ta + tb));
			if (insideBox(mid, box)) {
				out.add(new SegPiece(new MwvdGeom.Segment(a, b), a, b, edge.s1(), edge.s2()));
			}
		}
		return out;
	}

	private List<EdgePiece> clipArc(ArcPiece edge, Envelope box) {
		List<Coordinate> cuts = new ArrayList<>();
		addUnique(cuts, edge.start());
		addUnique(cuts, edge.end());

		for (MwvdGeom.Segment side : boxSides(box)) {
			for (Coordinate p : g.segmentArcIntersections(side, edge.arc())) {
				if (g.pointOnArc(p, edge.arc())) {
					addUnique(cuts, p);
				}
			}
		}

		cuts.sort(Comparator.comparingDouble(p -> arcParam(edge.arc(), p)));
		List<EdgePiece> out = new ArrayList<>();
		for (int i = 0; i + 1 < cuts.size(); i++) {
			Coordinate a = cuts.get(i);
			Coordinate b = cuts.get(i + 1);
			if (g.samePoint(a, b)) {
				continue;
			}
			double a0 = g.angle(edge.arc().circle().c(), a);
			double a1 = g.angle(edge.arc().circle().c(), b);
			double am = g.normAngle(a0 + 0.5 * g.ccwSweep(a0, a1));
			Coordinate mid = g.pointOnCircle(edge.arc().circle(), am);
			if (insideBox(mid, box)) {
				out.add(new ArcPiece(new MwvdGeom.Arc(edge.arc().circle(), a0, a1), a, b, edge.s1(), edge.s2()));
			}
		}
		return out;
	}

	private void addBoundaryEdges(HalfEdgeGraph graph, List<Coordinate> vertices, Envelope box) {
		addBoundarySide(graph, vertices, box, 0);
		addBoundarySide(graph, vertices, box, 1);
		addBoundarySide(graph, vertices, box, 2);
		addBoundarySide(graph, vertices, box, 3);
	}

	private void addBoundarySide(HalfEdgeGraph graph, List<Coordinate> vertices, Envelope box, int side) {
		List<Coordinate> sideVertices = new ArrayList<>();
		for (Coordinate p : vertices) {
			if (onBoundarySide(p, box, side)) {
				addUnique(sideVertices, p);
			}
		}
		sideVertices.sort(boundaryComparator(side));
		for (int i = 0; i + 1 < sideVertices.size(); i++) {
			Coordinate a = sideVertices.get(i);
			Coordinate b = sideVertices.get(i + 1);
			if (!g.samePoint(a, b)) {
				graph.addBoundary(a, b);
			}
		}
	}

	private Comparator<Coordinate> boundaryComparator(int side) {
		return switch (side) {
			case 0 -> Comparator.comparingDouble((Coordinate p) -> p.x).thenComparingDouble(p -> p.y);
			case 1 -> Comparator.comparingDouble((Coordinate p) -> p.y).thenComparingDouble(p -> p.x);
			case 2 -> Comparator.comparingDouble((Coordinate p) -> -p.x).thenComparingDouble(p -> p.y);
			case 3 -> Comparator.comparingDouble((Coordinate p) -> -p.y).thenComparingDouble(p -> p.x);
			default -> throw new IllegalArgumentException("side");
		};
	}

	private boolean onBoundarySide(Coordinate p, Envelope box, int side) {
		return switch (side) {
			case 0 -> g.eq(p.y, box.getMinY()) && g.le(box.getMinX(), p.x) && g.le(p.x, box.getMaxX());
			case 1 -> g.eq(p.x, box.getMaxX()) && g.le(box.getMinY(), p.y) && g.le(p.y, box.getMaxY());
			case 2 -> g.eq(p.y, box.getMaxY()) && g.le(box.getMinX(), p.x) && g.le(p.x, box.getMaxX());
			case 3 -> g.eq(p.x, box.getMinX()) && g.le(box.getMinY(), p.y) && g.le(p.y, box.getMaxY());
			default -> false;
		};
	}

	private List<MwvdGeom.Segment> boxSides(Envelope box) {
		Coordinate[] corners = boxCorners(box);
		return List.of(
				new MwvdGeom.Segment(corners[0], corners[1]),
				new MwvdGeom.Segment(corners[1], corners[2]),
				new MwvdGeom.Segment(corners[2], corners[3]),
				new MwvdGeom.Segment(corners[3], corners[0]));
	}

	private Coordinate[] boxCorners(Envelope box) {
		return new Coordinate[] {
				new Coordinate(box.getMinX(), box.getMinY()),
				new Coordinate(box.getMaxX(), box.getMinY()),
				new Coordinate(box.getMaxX(), box.getMaxY()),
				new Coordinate(box.getMinX(), box.getMaxY())
		};
	}

	private boolean insideBox(Coordinate p, Envelope box) {
		return g.le(box.getMinX(), p.x) && g.le(p.x, box.getMaxX())
				&& g.le(box.getMinY(), p.y) && g.le(p.y, box.getMaxY());
	}

	private boolean onBoundary(Coordinate p, Envelope box) {
		if (!insideBox(p, box)) {
			return false;
		}
		return g.eq(p.x, box.getMinX()) || g.eq(p.x, box.getMaxX())
				|| g.eq(p.y, box.getMinY()) || g.eq(p.y, box.getMaxY());
	}

	private double segmentParam(MwvdGeom.Segment seg, Coordinate p) {
		double dx = seg.b().x - seg.a().x;
		double dy = seg.b().y - seg.a().y;
		double denom = dx * dx + dy * dy;
		if (g.eq(denom, 0.0)) {
			return 0.0;
		}
		return ((p.x - seg.a().x) * dx + (p.y - seg.a().y) * dy) / denom;
	}

	private Coordinate segmentPoint(MwvdGeom.Segment seg, double t) {
		return new Coordinate(
				seg.a().x + (seg.b().x - seg.a().x) * t,
				seg.a().y + (seg.b().y - seg.a().y) * t);
	}

	private double arcParam(MwvdGeom.Arc arc, Coordinate p) {
		return g.ccwSweep(arc.a0(), g.angle(arc.circle().c(), p));
	}

	private void addUnique(List<Coordinate> out, Coordinate p) {
		for (Coordinate q : out) {
			if (g.samePoint(q, p)) {
				return;
			}
		}
		out.add(new Coordinate(p.x, p.y));
	}

	private EdgeShape toShape(EdgePiece edge) {
		return switch (edge) {
			case ArcPiece a -> {
				Coordinate[] forward = approxArc(a);
				yield new EdgeShape(
						a.start(),
						a.end(),
						forward,
						reversed(forward),
						g.normAngle(g.angle(a.arc().circle().c(), a.start()) + Math.PI / 2.0),
						g.normAngle(g.angle(a.arc().circle().c(), a.end()) - Math.PI / 2.0));
			}
			case SegPiece s -> new EdgeShape(
					s.start(),
					s.end(),
					new Coordinate[] { new Coordinate(s.start().x, s.start().y), new Coordinate(s.end().x, s.end().y) },
					new Coordinate[] { new Coordinate(s.end().x, s.end().y), new Coordinate(s.start().x, s.start().y) },
					g.angle(s.start(), s.end()),
					g.angle(s.end(), s.start()));
		};
	}

	private Coordinate[] reversed(Coordinate[] coords) {
		Coordinate[] out = new Coordinate[coords.length];
		for (int i = 0; i < coords.length; i++) {
			Coordinate c = coords[coords.length - 1 - i];
			out[i] = new Coordinate(c.x, c.y);
		}
		return out;
	}

	private Coordinate[] approxArc(ArcPiece edge) {
		double sweep = g.ccwSweep(edge.arc().a0(), edge.arc().a1());
		if (g.eq(sweep, 0.0)) {
			return new Coordinate[] {
					new Coordinate(edge.start().x, edge.start().y),
					new Coordinate(edge.end().x, edge.end().y)
			};
		}

		double r = edge.arc().circle().r();
		double maxStep = r <= chordTol ? Math.PI / 8.0
				: Math.min(Math.PI / 8.0, 2.0 * Math.acos(Math.max(-1.0, 1.0 - chordTol / r)));

		int n = Math.max(1, (int) Math.ceil(sweep / maxStep));
		Coordinate[] cs = new Coordinate[n + 1];
		cs[0] = new Coordinate(edge.start().x, edge.start().y);
		cs[n] = new Coordinate(edge.end().x, edge.end().y);

		for (int i = 1; i < n; i++) {
			double a = g.normAngle(edge.arc().a0() + sweep * i / n);
			cs[i] = g.pointOnCircle(edge.arc().circle(), a);
		}
		return cs;
	}

	private double signedArea(List<Coordinate> ring) {
		double area = 0.0;
		for (int i = 0; i + 1 < ring.size(); i++) {
			Coordinate a = ring.get(i);
			Coordinate b = ring.get(i + 1);
			area += a.x * b.y - b.x * a.y;
		}
		return 0.5 * area;
	}

	private void appendPath(List<Coordinate> ring, Coordinate[] coords) {
		for (Coordinate c : coords) {
			if (!ring.isEmpty() && g.samePoint(ring.get(ring.size() - 1), c)) {
				continue;
			}
			ring.add(new Coordinate(c.x, c.y));
		}
	}

	private record EdgeShape(
			Coordinate start,
			Coordinate end,
			Coordinate[] forward,
			Coordinate[] reverse,
			double forwardAngle,
			double reverseAngle) {
	}

	private final class HalfEdgeGraph {
		private final List<Vertex> vertices = new ArrayList<>();
		private final List<HalfEdge> halfEdges = new ArrayList<>();

		void add(EdgePiece edge) {
			EdgeShape shape = toShape(edge);
			add(shape.start(), shape.end(), shape.forward(), shape.reverse(), shape.forwardAngle(), shape.reverseAngle());
		}

		void addBoundary(Coordinate a, Coordinate b) {
			Coordinate[] forward = {
					new Coordinate(a.x, a.y),
					new Coordinate(b.x, b.y)
			};
			Coordinate[] reverse = {
					new Coordinate(b.x, b.y),
					new Coordinate(a.x, a.y)
			};
			add(a, b, forward, reverse, g.angle(a, b), g.angle(b, a));
		}

		private void add(Coordinate a, Coordinate b, Coordinate[] forward, Coordinate[] reverse, double forwardAngle, double reverseAngle) {
			if (g.samePoint(a, b)) {
				return;
			}
			Vertex from = vertex(a);
			Vertex to = vertex(b);

			HalfEdge f = new HalfEdge(from, to, forward, forwardAngle);
			HalfEdge r = new HalfEdge(to, from, reverse, reverseAngle);
			f.twin = r;
			r.twin = f;

			from.outgoing.add(f);
			to.outgoing.add(r);
			halfEdges.add(f);
			halfEdges.add(r);
		}

		private Vertex vertex(Coordinate p) {
			for (Vertex v : vertices) {
				if (g.samePoint(v.p, p)) {
					return v;
				}
			}
			Vertex v = new Vertex(new Coordinate(p.x, p.y));
			vertices.add(v);
			return v;
		}

		List<Polygon> walkFaces() {
			for (Vertex v : vertices) {
				v.outgoing.sort(Comparator.comparingDouble(e -> e.angle));
				for (int i = 0; i < v.outgoing.size(); i++) {
					v.outgoing.get(i).order = i;
				}
			}

			for (HalfEdge edge : halfEdges) {
				List<HalfEdge> out = edge.to.outgoing;
				if (out.isEmpty()) {
					continue;
				}
				int idx = edge.twin.order - 1;
				if (idx < 0) {
					idx += out.size();
				}
				edge.next = out.get(idx);
			}

			List<Polygon> faces = new ArrayList<>();
			for (HalfEdge start : halfEdges) {
				if (start.used) {
					continue;
				}

				List<Coordinate> ring = new ArrayList<>();
				HalfEdge edge = start;
				boolean closed = false;
				while (edge != null && !edge.used) {
					edge.used = true;
					appendPath(ring, edge.coords);
					edge = edge.next;
					if (edge == start) {
						closed = true;
						break;
					}
				}

				if (!closed || ring.size() < 3) {
					continue;
				}
				if (!g.samePoint(ring.get(0), ring.get(ring.size() - 1))) {
					ring.add(new Coordinate(ring.get(0).x, ring.get(0).y));
				}
				if (ring.size() < 4 || !g.lt(0.0, signedArea(ring))) {
					continue;
				}

				LinearRing shell = gf.createLinearRing(ring.toArray(Coordinate[]::new));
				if (!shell.isValid()) {
					continue;
				}
				Polygon face = gf.createPolygon(shell);
				if (!face.isEmpty() && face.isValid()) {
					faces.add(face);
				}
			}
			return faces;
		}
	}

	private static final class Vertex {
		private final Coordinate p;
		private final List<HalfEdge> outgoing = new ArrayList<>();

		private Vertex(Coordinate p) {
			this.p = p;
		}
	}

	private static final class HalfEdge {
		private final Vertex from;
		private final Vertex to;
		private final Coordinate[] coords;
		private final double angle;
		private HalfEdge twin;
		private HalfEdge next;
		private boolean used;
		private int order;

		private HalfEdge(Vertex from, Vertex to, Coordinate[] coords, double angle) {
			this.from = from;
			this.to = to;
			this.coords = coords;
			this.angle = angle;
		}
	}
}
