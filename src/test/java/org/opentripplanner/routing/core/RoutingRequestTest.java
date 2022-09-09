package org.opentripplanner.routing.core;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.opentripplanner.model.GenericLocation;
import org.opentripplanner.routing.api.request.RouteRequest;

public class RoutingRequestTest {

  @Test
  public void testRequest() {
    // TODO VIA: looks like some parts of this test are obsolete since method no longer exist

    RouteRequest request = new RouteRequest();
    //
    //    request.addMode(CAR);
    //    assertTrue(request.streetSubRequestModes.getCar());
    //    request.removeMode(CAR);
    //    assertFalse(request.streetSubRequestModes.getCar());

    request.setStreetSubRequestModes(new TraverseModeSet(TraverseMode.BICYCLE, TraverseMode.WALK));
    assertFalse(request.streetSubRequestModes.getCar());
    assertTrue(request.streetSubRequestModes.getBicycle());
    assertTrue(request.streetSubRequestModes.getWalk());
  }

  @Test
  public void testIntermediatePlaces() {
    // TODO VIA: those methods no longer exist (we will refactor them later). What should we do with this test?
    //    RoutingRequest req = new RoutingRequest();
    //    assertFalse(req.hasIntermediatePlaces());
    //
    //    req.clearIntermediatePlaces();
    //    assertFalse(req.hasIntermediatePlaces());
    //
    //    req.addIntermediatePlace(randomLocation());
    //    assertTrue(req.hasIntermediatePlaces());
    //
    //    req.clearIntermediatePlaces();
    //    assertFalse(req.hasIntermediatePlaces());
    //
    //    req.addIntermediatePlace(randomLocation());
    //    req.addIntermediatePlace(randomLocation());
    //    assertTrue(req.hasIntermediatePlaces());
  }

  @Test
  public void shouldCloneObjectFields() {
    var req = new RouteRequest();

    var clone = req.clone();

    assertNotSame(clone, req);
    assertNotSame(clone.raptorDebugging, req.raptorDebugging);

    assertEquals(50, req.numItineraries());
    assertEquals(50, clone.numItineraries());
  }

  private GenericLocation randomLocation() {
    return new GenericLocation(Math.random(), Math.random());
  }
}
