package org.opentripplanner.transit.model.site;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.opentripplanner.framework.i18n.I18NString;
import org.opentripplanner.framework.i18n.NonLocalizedString;
import org.opentripplanner.transit.model._data.TransitModelForTest;
import org.opentripplanner.transit.service.StopModel;

class GroupStopTest {

  private static TransitModelForTest TEST_MODEL = TransitModelForTest.of();

  private static final String ID = "1";
  private static final I18NString NAME = new NonLocalizedString("name");

  private static final StopLocation STOP_LOCATION = TEST_MODEL.stop("1:stop", 1d, 1d).build();
  private static final GroupStop subject = StopModel
    .of()
    .groupStop(TransitModelForTest.id(ID))
    .withName(NAME)
    .addLocation(STOP_LOCATION)
    .build();

  @Test
  void copy() {
    assertEquals(ID, subject.getId().getId());

    // Make a copy, and set the same name (nothing is changed)
    var copy = subject.copy().withName(NAME).build();

    assertSame(subject, copy);

    // Copy and change name
    copy = subject.copy().withName(new NonLocalizedString("v2")).build();

    // The two objects are not the same instance, but are equal(same id)
    assertNotSame(subject, copy);
    assertEquals(subject, copy);

    assertEquals(ID, copy.getId().getId());
    assertEquals(STOP_LOCATION, copy.getLocations().iterator().next());
    assertEquals("v2", copy.getName().toString());
  }

  @Test
  void sameAs() {
    assertTrue(subject.sameAs(subject.copy().build()));
    assertFalse(subject.sameAs(subject.copy().withId(TransitModelForTest.id("X")).build()));
    assertFalse(subject.sameAs(subject.copy().withName(new NonLocalizedString("X")).build()));
    assertFalse(
      subject.sameAs(
        subject.copy().addLocation(TransitModelForTest.of().stop("2:stop", 1d, 2d).build()).build()
      )
    );
  }
}
